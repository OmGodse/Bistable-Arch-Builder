"use strict";

const COLOR_STABLE = "#9C27B0";
const COLOR_UNSTABLE = "#FF0000"
const COLOR_FABRICATED = "#1976D2";

const HORIZONTAL_PADDING = 36;
const VERTICAL_PADDING = 24;

const WIDTH = 720 + 2 * HORIZONTAL_PADDING;
const HEIGHT = 360 + 2 * VERTICAL_PADDING;

let mode = "Draw";
let boundary = "Pinned";
let gridSpacing = 30;
const spacings = [20, 30, 36, 45, 60];

const drawTab = document.getElementById("tab-draw");
const codeTab = document.getElementById("tab-code");
const drawPanel = document.getElementById("panel-draw");
const codePanel = document.getElementById("panel-code");

const boundaryRadioButtons = document.getElementsByName("boundary-condition");
const slider = document.getElementById("slider");
const qInput = document.getElementById("q");

const graph = document.getElementById("graph");
const canvas = document.getElementById("canvas");
const results = document.getElementById("results");
const code = document.getElementById("python-code");

const lines = [];
const expression = [];
let forces = null;
let result = null;
let startCoordinates = null;
let highlightCoordinates = null;

const numberOfCoefficients = 30;

const RATIO = 5;
canvas.width = WIDTH * RATIO;
canvas.height = HEIGHT * RATIO;
canvas.style.width = WIDTH.toString() + "px";
canvas.style.height = HEIGHT.toString() + "px";
graph.style.width = (WIDTH - HORIZONTAL_PADDING * 2).toString() + "px";

const editorDiv = document.getElementById("editor");
const codeMirror = CodeMirror(editorDiv, {
    value: "const a = 2;\nreturn a * x * (1 - x);",
    mode: "javascript"
});

const context = canvas.getContext("2d");
invalidate();
invalidateTabLayout();

const worker = new Worker('Worker.js');
worker.onmessage = function (e) {
    code.innerHTML = e.data;
};

function getQ() {
    const defaultQ = 5;
    const input = qInput.value;
    if (input.length == 0) {
        return defaultQ;
    }
    const number = Number(input);
    if (Number.isNaN(number)) {
        return defaultQ;
    }
    return number;
}


function inEditMode() {
    return result == null && mode == "Draw";
}


function drawGrid() {
    context.fillText("1", -24, 0);
    context.fillText("0", -24, 180);
    context.fillText("1", 720 + 24, 180);

    context.beginPath();
    context.rect(0, 0, 720, 360);
    context.stroke();

    context.beginPath();
    for (var x = gridSpacing; x < 720; x += gridSpacing) {
        context.moveTo(x, 0);
        context.lineTo(x, 360);
    }
    for (var y = gridSpacing; y < 360; y += gridSpacing) {
        context.moveTo(0, y);
        context.lineTo(720, y);
    }
    context.stroke();
}

function drawPinnedBoundary() {
    const side = 12.5 * Math.sqrt(3);
    const radius = 4;
    drawTriangle(0, 180, side);
    drawTriangle(720, 180, side);
    drawCircle(0, 180, radius);
    drawCircle(720, 180, radius);
}

function drawTriangle(x, y, side) {
    const height = side * cos(PI / 6);
    const width = side * sin(PI / 6);

    context.beginPath();

    context.moveTo(x, y);
    context.lineTo(x + width, y + height);
    context.lineTo(x - width, y + height);
    context.lineTo(x, y);

    context.moveTo(x - width * 1.7, y + height);
    context.lineTo(x + width * 1.5, y + height);

    const dash = 10;
    const deltaX = dash * sin(PI / 6);
    const deltaY = dash * cos(PI / 6);

    const startPoints = [-width, -width / 2, 0, width / 2, width];
    startPoints.forEach((point) => {
        context.moveTo(x + point, y + height);
        context.lineTo(x + point - deltaX, y + height + deltaY);
    })

    context.stroke();
}

function drawCircle(x, y, radius) {
    context.beginPath();
    context.arc(x, y, radius, 0, 2 * Math.PI);
    context.stroke();

    context.beginPath();
    context.arc(x, y, radius - 1, 0, 2 * Math.PI);
    context.fill();
}

function drawFixedBoundary() {
    const height = 40;
    const dash = 15;
    const deltaX = dash * Math.sin(Math.PI / 3);
    const deltaY = dash * Math.cos(Math.PI / 3);

    context.beginPath();

    context.moveTo(0, 180 - height / 2);
    context.lineTo(0, 180 + height / 2);

    context.moveTo(720, 180 - height / 2);
    context.lineTo(720, 180 + height / 2);

    const startPoints = [180 - height / 2, 180 - height / 4, 180, 180 + height / 4, 180 + height / 2];
    startPoints.forEach((point) => {
        context.moveTo(0, point);
        context.lineTo(0 - deltaX, point + deltaY);

        context.moveTo(720, point);
        context.lineTo(720 + deltaX, point + deltaY);
    })

    context.stroke();
}


function plot(displayCoordinates) {
    context.beginPath();
    displayCoordinates.forEach((coordinate) => {
        const x = coordinate[0] * 720;
        const y = 180 - coordinate[1] * 180;
        if (x == 0) {
            context.moveTo(x, y);
        } else context.lineTo(x, y);
    });
    context.stroke();
}


function snap(offsetX, offsetY) {
    const X_LEFT = HORIZONTAL_PADDING;
    const X_RIGHT = WIDTH - HORIZONTAL_PADDING;
    const Y_TOP = VERTICAL_PADDING;
    const Y_BOTTOM = HEIGHT - VERTICAL_PADDING;

    let x;
    if (inRange(offsetX, X_LEFT, X_RIGHT)) {
        x = Math.round((offsetX - X_LEFT) / gridSpacing) * gridSpacing;
    }
    else {
        const xDiff = offsetX - X_LEFT;
        if (xDiff < 0) {
            x = 0;
        } else x = 720;
    }

    let y;
    if (inRange(offsetY, Y_TOP, Y_BOTTOM)) {
        y = Math.round((offsetY - Y_TOP) / gridSpacing) * gridSpacing;
    }
    else {
        const yDiff = offsetY - Y_TOP;
        if (yDiff < 0) {
            y = 0;
        } else y = 360;
    }

    return [x, y];
}



function reset() {
    lines.length = 0;
    expression.length = 0;
    result = null;
    startCoordinates = null;
    highlightCoordinates = null;
}

function invalidate() {
    context.reset();
    context.scale(RATIO, RATIO);
    context.clearRect(0, 0, WIDTH, HEIGHT);
    context.translate(HORIZONTAL_PADDING, VERTICAL_PADDING);

    context.lineWidth = 1;
    context.strokeStyle = "#BDBDBD";
    context.font = "Bold 16px Roboto";
    context.textAlign = "center";
    context.textBaseline = "middle";
    context.fillStyle = "#000000";
    drawGrid();

    context.lineWidth = 2;
    context.lineCap = "round";
    context.fillStyle = "#FFFFFF";
    context.strokeStyle = "#000000";
    if (boundary === "Pinned") {
        drawPinnedBoundary();
    }
    else {
        drawFixedBoundary();
    }

    if (result != null) {
        context.lineWidth = 2;
        context.strokeStyle = COLOR_FABRICATED;
        plot(result.fabricatedProfile.coordinates);

        result.toggledProfiles.forEach(profile => {
            if (profile.stable) {
                context.strokeStyle = COLOR_STABLE;
            }
            else {
                context.strokeStyle = COLOR_UNSTABLE;
                context.setLineDash([5, 5]);
            }
            plot(profile.coordinates);
        });

        if (result.forces != null) {
            const trace = {
                x: result.forces.x,
                y: result.forces.y,
                text: ["Switching Force (Asymmetric)", "Switching Force (Symmetric)", "Switchback Force (Symmetric)", "Switchback Force (Asymmetric)"],
                mode: "markers",
                type: "scatter",
                marker: {
                    size: 8
                },
                hoverlabel: {
                    font: { family: "Roboto", size: 14 }
                }
            };
            const layout = {
                font: {
                    family: "Roboto",
                    size: 14
                },
                xaxis: {
                    title: {
                        text: "<b><span style=\"letter-spacing:0.5px\">Displacement</span><b>",
                        font: { size: 16, color: "black" }
                    },
                    rangemode: "tozero",
                    fixedrange: true
                },
                yaxis: {
                    title: {
                        text: "<b><span style=\"letter-spacing:0.5px\">Force</span><b>",
                        font: { size: 16, color: "black" }
                    },
                    fixedrange: true
                },
                margin: { t: 16, r: 0, pad: 8 }
            };
            const options = { displaylogo: false };
            Plotly.newPlot(graph, [trace], layout, options);
            graph.style.display = "block";
        }
        else {
            graph.style.display = "none";
        }

        worker.postMessage(result);
        results.style.display = "block";
    }
    else {
        results.style.display = "none";
        graph.style.display = "none";

        if (mode === "Draw") {
            lines.forEach((line) => {
                context.beginPath();
                context.moveTo(line[0], line[1]);
                context.lineTo(line[2], line[3]);
                context.stroke();
            });

            if (startCoordinates != null) {
                context.strokeStyle = "#FF0000";
                context.setLineDash([5, 5]);

                context.beginPath();
                context.moveTo(startCoordinates[0], startCoordinates[1]);
                context.lineTo(highlightCoordinates[0], highlightCoordinates[1]);
                context.stroke();
            }

            if (highlightCoordinates != null) {
                context.beginPath();
                context.arc(highlightCoordinates[0], highlightCoordinates[1], 5, 0, 2 * Math.PI);
                context.fillStyle = "#FF0000";
                context.fill();
            }
        }
        else {
            plot(expression);
        }
    }

    if (inEditMode()) {
        canvas.style.cursor = "crosshair";
    }
    else {
        canvas.style.cursor = "auto";
    }
}

function invalidateTabLayout() {
    if (mode === "Draw") {
        drawTab.className = "selected-tab";
        drawPanel.style.display = "block";
        codeTab.className = "tab";
        codePanel.style.display = "none";
    }
    else {
        codeTab.className = "selected-tab";
        codePanel.style.display = "block";
        drawTab.className = "tab";
        drawPanel.style.display = "none";
    }
}


function addExpression(string) {
    const functionA = new Function("x", string);
    const SCALE = 1000;
    for (let scaledX = 0; scaledX <= SCALE; scaledX++) {
        const x = scaledX / SCALE;
        const y = functionA(x);
        expression[scaledX] = [x, y];
    }

    const Q = getQ();
    const initialCoefficients = calculateCoefficientsExpression(numberOfCoefficients, expression, boundary);
    const hMid = calculateHMidExpression(expression);
    determineProfileType(initialCoefficients, Q, hMid, boundary);
}

function calculateFourierSeries() {
    if (continuous(lines)) {
        const Q = getQ();
        const initialCoefficients = calculateCoefficientsLines(numberOfCoefficients, lines, boundary)
        const hMid = calculateHMid(lines);
        determineProfileType(initialCoefficients, Q, hMid, boundary);
    }
}


function determineProfileType(initialCoefficients, Q, hMid, boundary) {
    result = null;

    if (hMid > 0) {
        calculateToggledProfiles(initialCoefficients, Q, boundary);
    }
    else {
        calculateFabricatedProfile(initialCoefficients, Q, boundary);
    }
}

function calculateToggledProfiles(hCoefficients, Q, boundary) {
    const fabricatedCoordinates = calculateCoordinates(hCoefficients, boundary);
    const fabricatedProfile = { coefficients: hCoefficients, coordinates: fabricatedCoordinates }

    const A1s = calculateA1s(hCoefficients, Q, boundary);

    const toggledProfiles = [];

    let stableWCoefficients = null;
    A1s.forEach(A1 => {
        const toggledCoefficients = calculateWCoefficients(hCoefficients, A1, boundary);
        const toggledCoordinates = calculateCoordinates(toggledCoefficients, boundary);
        const hessian = calculateHessian(hCoefficients, toggledCoefficients, Q, boundary);
        const eigenvalues = math.eigs(hessian, { eigenvectors: false }).values;
        const stable = eigenvalues.every((eigenvalue) => eigenvalue >= 0);
        if (stable) {
            stableWCoefficients = toggledCoefficients;
        }
        toggledProfiles.push({ coefficients: toggledCoefficients, coordinates: toggledCoordinates, stable: stable })
    })

    let forces = null;
    if (stableWCoefficients != null) {
        forces = calculateForces(hCoefficients, stableWCoefficients, Q, boundary);
    }

    result = { fabricatedProfile: fabricatedProfile, toggledProfiles: toggledProfiles, forces: forces };
}

function calculateFabricatedProfile(wCoefficients, Q, boundary) {
    const toggledCoordinates = calculateCoordinates(wCoefficients, boundary);
    const hCoefficients = calculateHCoefficients(wCoefficients, Q, boundary);
    const fabricatedCoordinates = calculateCoordinates(hCoefficients, boundary);
    const hessian = calculateHessian(hCoefficients, wCoefficients, Q, boundary);

    const eigenvalues = math.eigs(hessian, { eigenvectors: false }).values
    const stable = eigenvalues.every((eigenvalue) => eigenvalue >= 0);
    let forces = null;
    if (stable) {
        forces = calculateForces(hCoefficients, wCoefficients, Q, boundary);
    }

    const toggledProfiles = [{ coefficients: wCoefficients, coordinates: toggledCoordinates, stable: stable }]
    const fabricatedProfile = { coefficients: hCoefficients, coordinates: fabricatedCoordinates }

    result = { fabricatedProfile: fabricatedProfile, toggledProfiles: toggledProfiles, forces: forces };
}


drawTab.onclick = function () {
    mode = "Draw";
    invalidateTabLayout();
    invalidate();
}

codeTab.onclick = function () {
    mode = "Code";
    invalidateTabLayout();
    invalidate();
}

boundaryRadioButtons.forEach((button) => {
    button.onchange = function () {
        boundary = button.value;
        invalidate();
    }
});

slider.oninput = function () {
    gridSpacing = spacings[slider.value];
    invalidate();
}

document.getElementById("reset").onclick = function () {
    reset();
    invalidate();
};

document.getElementById("calculate").onclick = function () {
    addExpression(codeMirror.getValue());
    invalidate();
}

canvas.oncontextmenu = function (event) {
    if (startCoordinates != null) {
        event.preventDefault();
        startCoordinates = null;
        invalidate();
    }
};

canvas.onclick = function (event) {
    if (inEditMode()) {
        const coordinates = snap(event.offsetX, event.offsetY);
        if (startCoordinates != null) {
            let x1;
            let y1;
            let x2;
            let y2;
            if (startCoordinates[0] < coordinates[0]) {
                x1 = startCoordinates[0];
                y1 = startCoordinates[1];
                x2 = coordinates[0];
                y2 = coordinates[1];
            }
            else {
                x1 = coordinates[0];
                y1 = coordinates[1];
                x2 = startCoordinates[0];
                y2 = startCoordinates[1];
            }
            lines.push([x1, y1, x2, y2]);
            startCoordinates = null;
            calculateFourierSeries();
        }
        else {
            startCoordinates = coordinates;
        }
        invalidate();
    }
};

canvas.onmousemove = function (event) {
    if (inEditMode()) {
        highlightCoordinates = snap(event.offsetX, event.offsetY);
        invalidate();
    }
};

canvas.onmouseleave = function () {
    if (inEditMode()) {
        if (startCoordinates == null) {
            highlightCoordinates = null;
            invalidate();
        }
    }
};