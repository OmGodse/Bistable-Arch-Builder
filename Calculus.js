"use strict";

const sin = Math.sin;
const cos = Math.cos;
const pow = Math.pow;
const abs = Math.abs;
const PI = Math.PI;

function continuous(lines) {
    let xStart = 0;
    let yStart = 180;
    let numberOfConnectingLines = 0;
    while (true) {
        let indices = findIndicesForLines(xStart, yStart, lines);
        if (indices.length < 1) {
            console.log("Gap at (" + xStart + ", " + yStart + ")");
            break;
        }
        else if (indices.length > 1) {
            console.log("Ambiguity at (" + xStart + ", " + yStart + ")");
            break;
        }
        else {
            const index = indices[0];
            const line = lines[index];
            xStart = line[2];
            yStart = line[3];
            numberOfConnectingLines++;
        }
    }

    if (numberOfConnectingLines < lines.length) {
        console.log("Extra lines");
        return false;
    }
    else return xStart === 720 && yStart === 180;
}

function findIndicesForLines(xStart, yStart, lines) {
    const indices = [];
    lines.forEach((line, index) => {
        const x1 = line[0];
        const y1 = line[1];
        if (x1 === xStart && y1 === yStart) {
            indices.push(index);
        }
    });
    return indices;
}


function calculateHMid(lines) {
    let hMid = 0;
    lines.forEach((line) => {
        const y1 = 1 - line[1] / 180;
        const y2 = 1 - line[3] / 180;
        if (abs(y1) > abs(hMid)) {
            hMid = y1;
        }
        if (abs(y2) > abs(hMid)) {
            hMid = y2;
        }
    });
    return hMid;
}

function calculateHMidExpression(coordinates) {
    let hMid = 0;
    coordinates.forEach((coordinate) => {
        const y = coordinate[1];
        if (abs(y) > abs(hMid)) {
            hMid = y;
        }
    })
    return hMid;
}


function fourierCoefficientIntegralLines(i, lines, antiderivative) {
    let integral = 0;
    lines.forEach((line) => {
        if (line[0] != line[2]) {
            const x1 = line[0] / 720;
            const x2 = line[2] / 720;
            const y1 = 1 - line[1] / 180;
            const y2 = 1 - line[3] / 180;
            const m = (y2 - y1) / (x2 - x1);
            const c = - (m * x1) + y1;
            integral += antiderivative(x2, m, c, i) - antiderivative(x1, m, c, i);
        }
    })
    return integral;
}

function pinned_lineAntiderivative(x, m, c, i) {
    return (m * sin(PI * i * x) - PI * i * (m * x + c) * cos(PI * i * x)) / (PI * PI * i * i);
}

function fixed_lineAntiderivative(x, m, c, i) {
    const N = M(i, "Fixed");
    if (i % 2 == 0) {
        return - (m * (N * N * x - 2) + N * N * c) * sin(N * x) / (N * N * N) - (m * (2 * x + 1) + 2 * c) * cos(N * x) / (N * N) - 2 * m * x * x * x / 3 + (m - 2 * c) * x * x / 2 + c * x;
    }
    else return - (m * x + c) * sin(N * x) / N - m * cos(N * x) / (N * N) + m * x * x / 2 + c * x
}


function calculateCoefficientsLines(n, lines, boundary) {
    if (boundary === "Pinned") {
        return pinned_calculateCoefficients(n, (i) => fourierCoefficientIntegralLines(i, lines, pinned_lineAntiderivative));
    }
    else {
        return fixed_calculateCoefficients(n, (i) => fourierCoefficientIntegralLines(i, lines, fixed_lineAntiderivative));
    }
}

function calculateCoefficientsExpression(n, coordinates, boundary) {
    if (boundary === "Pinned") {
        return pinned_calculateCoefficients(n, (i) => {
            const basis = (x) => sin(PI * i * x);
            return numericallyIntegrate(coordinates, basis)
        });
    }
    else {
        return fixed_calculateCoefficients(n, (i) => {
            const m = M(i, "Fixed");
            let basis;
            if (i % 2 == 0) {
                basis = (x) => 1 - 2 * x - cos(m * x) + 2 * sin(m * x) / m;
            }
            else basis = (x) => 1 - cos(m * x);
            return numericallyIntegrate(coordinates, basis);
        });
    }
}


function pinned_calculateCoefficients(n, coefficientIntegral) {
    const coefficients = Array(n);
    for (let i = 1; i <= n; i++) {
        const integral = coefficientIntegral(i);
        if (abs(integral) < 1e-10) {
            coefficients[i - 1] = 0;
        }
        else {
            coefficients[i - 1] = 2 * integral;
        }
    }
    return coefficients;
}

function fixed_calculateCoefficients(n, coefficientIntegral) {
    const A = Array(n);
    const b = Array(n);
    for (let i = 1; i <= n; i++) {
        const row = Array(n);
        for (let j = 1; j <= n; j++) {
            let coefficient;
            if (i % 2 == 0) {
                if (j % 2 == 0) {
                    if (j == i) {
                        coefficient = 5 / 6;
                    }
                    else coefficient = 1 / 3;
                }
                else coefficient = 0;
            }
            else {
                if (j % 2 == 0) {
                    coefficient = 0;
                }
                else {
                    if (j == i) {
                        coefficient = 3 / 2;
                    }
                    else coefficient = 1;
                }
            }
            row[j - 1] = coefficient;
        }
        A[i - 1] = row;
        const integral = coefficientIntegral(i);
        if (abs(integral) < 1e-10) {
            b[i - 1] = 0;
        }
        else {
            b[i - 1] = integral;
        }
    }
    const solution = math.lusolve(A, b);
    return solution.map((array) => array[0]);
}


function numericallyIntegrate(coordinates, basis) {
    const n = coordinates.length - 1;
    const fA = coordinates[0][1] * basis(0);
    const fB = coordinates[n][1] * basis(1);
    let integral = 0.5 * (fA + fB);
    for (let i = 1; i < n; i++) {
        integral += coordinates[i][1] * basis(coordinates[i][0]);
    }
    return integral / n;
}


function M(i, boundary) {
    if (boundary === "Pinned") {
        return i * PI;
    }
    else {
        if (i % 2 == 0) {
            const index = (i / 2) - 1;
            const roots = [8.986818915818128, 15.450503673875414, 21.808243318857798, 28.132387825649133, 34.44151054384522, 40.74260591857664, 47.03890499737044, 53.332108517622466, 59.62319758173693, 65.91277807964462, 72.20124448875083, 78.48886472222652, 84.77582713620333, 91.06226802798213, 97.3482884639088];
            return roots[index];
        }
        else {
            return (i + 1) * PI;
        }
    }
}


function A1Equation(x, hCoefficients, Q, boundary) {
    let equationSeries = 0;
    hCoefficients.forEach((a, index) => {
        const i = index + 1;
        const c = pow(M(1, boundary) / M(i, boundary), 2);
        equationSeries += a * a * (x * c - 2) / pow(1 - x * c, 2);
    });
    return 3 * Q * Q * equationSeries - 1;
}

function calculateHCoefficients(wCoefficients, Q, boundary) {
    let firstSeries = 0;
    let secondSeries = 0;
    wCoefficients.forEach((coefficient, index) => {
        const i = index + 1;
        firstSeries += pow(M(1, boundary) * coefficient / M(i, boundary), 2)
        secondSeries += coefficient * coefficient;
    });

    const A1 = wCoefficients[0];
    const a1 = (A1 * (firstSeries - 2 * secondSeries - 1 / (3 * Q * Q))) / firstSeries;

    return wCoefficients.map((coefficient, index) => {
        if (index > 0) {
            const i = index + 1
            return coefficient * (1 - (1 - a1 / A1) * pow(M(1, boundary) / M(i, boundary), 2));
        }
        else return a1;
    });
}

function calculateWCoefficients(hCoefficients, A1, boundary) {
    const a1 = hCoefficients[0];
    return hCoefficients.map((coefficient, index) => {
        if (index > 0) {
            const i = index + 1;
            return coefficient / (1 - (1 - a1 / A1) * pow(M(1, boundary) / M(i, boundary), 2));
        }
        else return A1;
    });
}

function calculateHessian(hCoefficients, wCoefficients, Q, boundary) {
    const length = hCoefficients.length;

    let C = 0;
    for (let i = 1; i <= length; i++) {
        const a = hCoefficients[i - 1];
        const A = wCoefficients[i - 1];
        const Mi = M(i, boundary);
        C += Mi * Mi * (a * a - A * A);
    }

    const hessian = Array(length);
    for (let i = 1; i <= length; i++) {
        const row = Array(length);
        for (let j = 1; j <= length; j++) {
            const Ai = wCoefficients[i - 1];
            const Mi = M(i, boundary);
            let H;
            if (i != j) {
                const Aj = wCoefficients[j - 1];
                const Mj = M(j, boundary);
                H = 3 * Q * Q * Ai * Aj * Mi * Mi * Mj * Mj;
            }
            else {
                H = pow(Mi, 4) / 2 - (3 * Q * Q * Mi * Mi * C) / 2 + 3 * Q * Q * Ai * Ai * pow(Mi, 4);
            }
            row[j - 1] = H;
        }
        hessian[i - 1] = row;
    }

    return hessian;
}

function calculateA1s(hCoefficients, Q, boundary) {
    //console.time("New");
    const tolerance = 1e-6;
    const n = 200;
    const array = linspace(-10, 10, n + 1);
    console.log(array);
    const luckyRoots = [];
    const intervals = [];
    for (let i = 0; i < n; i++) {
        const a = array[i];
        const b = array[i + 1];
        const fA = A1Equation(a, hCoefficients, Q, boundary);
        const fB = A1Equation(b, hCoefficients, Q, boundary);
        if (abs(fA) < tolerance) {
            luckyRoots.push(a);
        }
        if (abs(fB) < tolerance) {
            luckyRoots.push(b);
        }
        if (fA * fB < 0) {
            intervals.push([a, b]);
        }
    }
    const roots = [];
    intervals.forEach((interval) => {
        const a = interval[0];
        const b = interval[1];
        //const fA = A1Equation(a, hCoefficients, Q);
        //const fB = A1Equation(b, hCoefficients, Q);
        //console.log(a, b, fA, fB);
        const root = bisect(hCoefficients, Q, boundary, a, b, 30, tolerance);
        if (root != undefined) {
            roots.push(root);
        }
    });
    //console.timeEnd("New");
    //console.log(luckyRoots);
    return roots.map((x) => hCoefficients[0] / (1 - x));
}


function calculateCoordinates(coefficients, boundary) {
    const SCALE = 1000;
    const coordinates = Array(SCALE + 1);
    for (let scaledX = 0; scaledX <= SCALE; scaledX++) {
        const x = scaledX / SCALE;
        const y = calculateYCoordinate(coefficients, boundary, x);
        coordinates[scaledX] = [x, y];
    }
    return coordinates;
}

function calculateYCoordinate(coefficients, boundary, x) {
    let y = 0;
    coefficients.forEach((coefficient, index) => {
        const i = index + 1;
        const m = M(i, boundary);
        if (boundary === "Pinned") {
            y += coefficient * sin(m * x)
        }
        else {
            if (i % 2 == 0) {
                y += coefficient * (1 - 2 * x - cos(m * x) + 2 * sin(m * x) / m);
            }
            else {
                y += coefficient * (1 - cos(m * x));
            }
        }
    });
    return y;
}


function FSymmetric(x, parameters) {
    const A1 = x[0];
    const A2 = x[1];
    const A3 = x[2];
    const lambda1 = x[3];
    const lambda2 = x[4];

    const Q = parameters[0];
    const a1 = parameters[1];
    const a2 = parameters[2];
    const a3 = parameters[3];

    const M1 = PI;
    const M2 = 2 * PI;
    const M3 = 3 * PI;

    const equation = [
        -3 * A1 * A2 * pow(M1, 2) * pow(M2, 2) * pow(Q, 2) * lambda2 + 3 * A1 * A3 * pow(M1, 2) * pow(M3, 2) * pow(Q, 2) - lambda1 * (3 * pow(A1, 2) * pow(M1, 4) * pow(Q, 2) + 3 * A1 * A3 * pow(M1, 2) * pow(M3, 2) * pow(Q, 2) + (1 / 2) * pow(M1, 4) - 3 / 2 * pow(M1, 2) * pow(Q, 2) * (-pow(A1, 2) * pow(M1, 2) - pow(A2, 2) * pow(M2, 2) - pow(A3, 2) * pow(M3, 2) + pow(M1, 2) * pow(a1, 2) + pow(M2, 2) * pow(a2, 2) + pow(M3, 2) * pow(a3, 2))),
        3 * A2 * A3 * pow(M2, 2) * pow(M3, 2) * pow(Q, 2) - lambda1 * (3 * A1 * A2 * pow(M1, 2) * pow(M2, 2) * pow(Q, 2) + 3 * A2 * A3 * pow(M2, 2) * pow(M3, 2) * pow(Q, 2)) - lambda2 * (3 * pow(A2, 2) * pow(M2, 4) * pow(Q, 2) + (1 / 2) * pow(M2, 4) - 3 / 2 * pow(M2, 2) * pow(Q, 2) * (-pow(A1, 2) * pow(M1, 2) - pow(A2, 2) * pow(M2, 2) - pow(A3, 2) * pow(M3, 2) + pow(M1, 2) * pow(a1, 2) + pow(M2, 2) * pow(a2, 2) + pow(M3, 2) * pow(a3, 2))),
        -3 * A2 * A3 * pow(M2, 2) * pow(M3, 2) * pow(Q, 2) * lambda2 + (1 / 2) * pow(M3, 2) * (3 * pow(A1, 2) * pow(M1, 2) * pow(Q, 2) + 3 * pow(A2, 2) * pow(M2, 2) * pow(Q, 2) + 9 * pow(A3, 2) * pow(M3, 2) * pow(Q, 2) - 3 * pow(M1, 2) * pow(Q, 2) * pow(a1, 2) - 3 * pow(M2, 2) * pow(Q, 2) * pow(a2, 2) - 3 * pow(M3, 2) * pow(Q, 2) * pow(a3, 2) + pow(M3, 2)) - lambda1 * (3 * A1 * A3 * pow(M1, 2) * pow(M3, 2) * pow(Q, 2) + (1 / 2) * pow(M3, 2) * (3 * pow(A1, 2) * pow(M1, 2) * pow(Q, 2) + 3 * pow(A2, 2) * pow(M2, 2) * pow(Q, 2) + 9 * pow(A3, 2) * pow(M3, 2) * pow(Q, 2) - 3 * pow(M1, 2) * pow(Q, 2) * pow(a1, 2) - 3 * pow(M2, 2) * pow(Q, 2) * pow(a2, 2) - 3 * pow(M3, 2) * pow(Q, 2) * pow(a3, 2) + pow(M3, 2))),
        (3 / 2) * A1 * pow(M1, 2) * pow(Q, 2) * (-pow(A1, 2) * pow(M1, 2) - pow(A2, 2) * pow(M2, 2) - pow(A3, 2) * pow(M3, 2) + pow(M1, 2) * pow(a1, 2) + pow(M2, 2) * pow(a2, 2) + pow(M3, 2) * pow(a3, 2)) - 1 / 2 * pow(M1, 4) * (A1 - a1) - 1 / 2 * pow(M3, 2) * (3 * pow(A1, 2) * A3 * pow(M1, 2) * pow(Q, 2) + 3 * pow(A2, 2) * A3 * pow(M2, 2) * pow(Q, 2) + 3 * pow(A3, 3) * pow(M3, 2) * pow(Q, 2) - 3 * A3 * pow(M1, 2) * pow(Q, 2) * pow(a1, 2) - 3 * A3 * pow(M2, 2) * pow(Q, 2) * pow(a2, 2) - 3 * A3 * pow(M3, 2) * pow(Q, 2) * pow(a3, 2) + A3 * pow(M3, 2) - pow(M3, 2) * a3),
        (3 / 2) * A2 * pow(M2, 2) * pow(Q, 2) * (-pow(A1, 2) * pow(M1, 2) - pow(A2, 2) * pow(M2, 2) - pow(A3, 2) * pow(M3, 2) + pow(M1, 2) * pow(a1, 2) + pow(M2, 2) * pow(a2, 2) + pow(M3, 2) * pow(a3, 2)) - 1 / 2 * pow(M2, 4) * (A2 - a2)
    ]
    const jacobian = [
        [
            -3 * A2 * pow(M1, 2) * pow(M2, 2) * pow(Q, 2) * lambda2 + 3 * A3 * pow(M1, 2) * pow(M3, 2) * pow(Q, 2) - lambda1 * (9 * A1 * pow(M1, 4) * pow(Q, 2) + 3 * A3 * pow(M1, 2) * pow(M3, 2) * pow(Q, 2)),
            -3 * A1 * pow(M1, 2) * pow(M2, 2) * pow(Q, 2) * lambda2 - 3 * A2 * pow(M1, 2) * pow(M2, 2) * pow(Q, 2) * lambda1,
            3 * A1 * pow(M1, 2) * pow(M3, 2) * pow(Q, 2) - lambda1 * (3 * A1 * pow(M1, 2) * pow(M3, 2) * pow(Q, 2) + 3 * A3 * pow(M1, 2) * pow(M3, 2) * pow(Q, 2)),
            -3 * pow(A1, 2) * pow(M1, 4) * pow(Q, 2) - 3 * A1 * A3 * pow(M1, 2) * pow(M3, 2) * pow(Q, 2) - 1 / 2 * pow(M1, 4) + (3 / 2) * pow(M1, 2) * pow(Q, 2) * (-pow(A1, 2) * pow(M1, 2) - pow(A2, 2) * pow(M2, 2) - pow(A3, 2) * pow(M3, 2) + pow(M1, 2) * pow(a1, 2) + pow(M2, 2) * pow(a2, 2) + pow(M3, 2) * pow(a3, 2)),
            -3 * A1 * A2 * pow(M1, 2) * pow(M2, 2) * pow(Q, 2)
        ],
        [
            -3 * A1 * pow(M1, 2) * pow(M2, 2) * pow(Q, 2) * lambda2 - 3 * A2 * pow(M1, 2) * pow(M2, 2) * pow(Q, 2) * lambda1,
            -9 * A2 * pow(M2, 4) * pow(Q, 2) * lambda2 + 3 * A3 * pow(M2, 2) * pow(M3, 2) * pow(Q, 2) - lambda1 * (3 * A1 * pow(M1, 2) * pow(M2, 2) * pow(Q, 2) + 3 * A3 * pow(M2, 2) * pow(M3, 2) * pow(Q, 2)),
            -3 * A2 * pow(M2, 2) * pow(M3, 2) * pow(Q, 2) * lambda1 + 3 * A2 * pow(M2, 2) * pow(M3, 2) * pow(Q, 2) - 3 * A3 * pow(M2, 2) * pow(M3, 2) * pow(Q, 2) * lambda2,
            -3 * A1 * A2 * pow(M1, 2) * pow(M2, 2) * pow(Q, 2) - 3 * A2 * A3 * pow(M2, 2) * pow(M3, 2) * pow(Q, 2),
            -3 * pow(A2, 2) * pow(M2, 4) * pow(Q, 2) - 1 / 2 * pow(M2, 4) + (3 / 2) * pow(M2, 2) * pow(Q, 2) * (-pow(A1, 2) * pow(M1, 2) - pow(A2, 2) * pow(M2, 2) - pow(A3, 2) * pow(M3, 2) + pow(M1, 2) * pow(a1, 2) + pow(M2, 2) * pow(a2, 2) + pow(M3, 2) * pow(a3, 2))
        ],
        [
            3 * A1 * pow(M1, 2) * pow(M3, 2) * pow(Q, 2) - lambda1 * (3 * A1 * pow(M1, 2) * pow(M3, 2) * pow(Q, 2) + 3 * A3 * pow(M1, 2) * pow(M3, 2) * pow(Q, 2)),
            -3 * A2 * pow(M2, 2) * pow(M3, 2) * pow(Q, 2) * lambda1 + 3 * A2 * pow(M2, 2) * pow(M3, 2) * pow(Q, 2) - 3 * A3 * pow(M2, 2) * pow(M3, 2) * pow(Q, 2) * lambda2,
            -3 * A2 * pow(M2, 2) * pow(M3, 2) * pow(Q, 2) * lambda2 + 9 * A3 * pow(M3, 4) * pow(Q, 2) - lambda1 * (3 * A1 * pow(M1, 2) * pow(M3, 2) * pow(Q, 2) + 9 * A3 * pow(M3, 4) * pow(Q, 2)),
            -3 * A1 * A3 * pow(M1, 2) * pow(M3, 2) * pow(Q, 2) - 1 / 2 * pow(M3, 2) * (3 * pow(A1, 2) * pow(M1, 2) * pow(Q, 2) + 3 * pow(A2, 2) * pow(M2, 2) * pow(Q, 2) + 9 * pow(A3, 2) * pow(M3, 2) * pow(Q, 2) - 3 * pow(M1, 2) * pow(Q, 2) * pow(a1, 2) - 3 * pow(M2, 2) * pow(Q, 2) * pow(a2, 2) - 3 * pow(M3, 2) * pow(Q, 2) * pow(a3, 2) + pow(M3, 2)),
            -3 * A2 * A3 * pow(M2, 2) * pow(M3, 2) * pow(Q, 2)
        ],
        [
            -3 * pow(A1, 2) * pow(M1, 4) * pow(Q, 2) - 3 * A1 * A3 * pow(M1, 2) * pow(M3, 2) * pow(Q, 2) - 1 / 2 * pow(M1, 4) + (3 / 2) * pow(M1, 2) * pow(Q, 2) * (-pow(A1, 2) * pow(M1, 2) - pow(A2, 2) * pow(M2, 2) - pow(A3, 2) * pow(M3, 2) + pow(M1, 2) * pow(a1, 2) + pow(M2, 2) * pow(a2, 2) + pow(M3, 2) * pow(a3, 2)),
            -3 * A1 * A2 * pow(M1, 2) * pow(M2, 2) * pow(Q, 2) - 3 * A2 * A3 * pow(M2, 2) * pow(M3, 2) * pow(Q, 2),
            -3 * A1 * A3 * pow(M1, 2) * pow(M3, 2) * pow(Q, 2) - 1 / 2 * pow(M3, 2) * (3 * pow(A1, 2) * pow(M1, 2) * pow(Q, 2) + 3 * pow(A2, 2) * pow(M2, 2) * pow(Q, 2) + 9 * pow(A3, 2) * pow(M3, 2) * pow(Q, 2) - 3 * pow(M1, 2) * pow(Q, 2) * pow(a1, 2) - 3 * pow(M2, 2) * pow(Q, 2) * pow(a2, 2) - 3 * pow(M3, 2) * pow(Q, 2) * pow(a3, 2) + pow(M3, 2)),
            0,
            0
        ],
        [
            -3 * A1 * A2 * pow(M1, 2) * pow(M2, 2) * pow(Q, 2),
            -3 * pow(A2, 2) * pow(M2, 4) * pow(Q, 2) - 1 / 2 * pow(M2, 4) + (3 / 2) * pow(M2, 2) * pow(Q, 2) * (-pow(A1, 2) * pow(M1, 2) - pow(A2, 2) * pow(M2, 2) - pow(A3, 2) * pow(M3, 2) + pow(M1, 2) * pow(a1, 2) + pow(M2, 2) * pow(a2, 2) + pow(M3, 2) * pow(a3, 2)),
            -3 * A2 * A3 * pow(M2, 2) * pow(M3, 2) * pow(Q, 2),
            0,
            0
        ]
    ];

    return [equation, jacobian];
}

function FAsymmetric(x, parameters) {
    const A1 = x[0];
    const A2 = x[1];
    const A3 = x[2];
    const F = x[3];

    const Q = parameters[0];
    const a1 = parameters[1];
    const a2 = parameters[2];
    const a3 = parameters[3];

    const M1 = PI;
    const M2 = 2 * PI;
    const M3 = 3 * PI;

    const equation = [
        -3 / 2 * A1 * pow(M1, 2) * pow(Q, 2) * (-pow(A1, 2) * pow(M1, 2) - pow(A2, 2) * pow(M2, 2) - pow(A3, 2) * pow(M3, 2) + pow(M1, 2) * pow(a1, 2) + pow(M2, 2) * pow(a2, 2) + pow(M3, 2) * pow(a3, 2)) + F + (1 / 2) * pow(M1, 4) * (A1 - a1),
        3 * A1 * A2 * pow(M1, 2) * pow(M2, 2) * pow(Q, 2),
        3 * pow(A2, 2) * pow(M2, 4) * pow(Q, 2) + (1 / 2) * pow(M2, 4) - 3 / 2 * pow(M2, 2) * pow(Q, 2) * (-pow(A1, 2) * pow(M1, 2) - pow(A2, 2) * pow(M2, 2) - pow(A3, 2) * pow(M3, 2) + pow(M1, 2) * pow(a1, 2) + pow(M2, 2) * pow(a2, 2) + pow(M3, 2) * pow(a3, 2)),
        3 * A2 * A3 * pow(M2, 2) * pow(M3, 2) * pow(Q, 2),
        -3 / 2 * A3 * pow(M3, 2) * pow(Q, 2) * (-pow(A1, 2) * pow(M1, 2) - pow(A2, 2) * pow(M2, 2) - pow(A3, 2) * pow(M3, 2) + pow(M1, 2) * pow(a1, 2) + pow(M2, 2) * pow(a2, 2) + pow(M3, 2) * pow(a3, 2)) - F + (1 / 2) * pow(M3, 4) * (A3 - a3)
    ]
    const row1 = [
        3 * Math.pow(A1, 2) * Math.pow(M1, 4) * Math.pow(Q, 2) + (1 / 2) * Math.pow(M1, 4) - 3 / 2 * Math.pow(M1, 2) * Math.pow(Q, 2) * (-Math.pow(A1, 2) * Math.pow(M1, 2) - Math.pow(A2, 2) * Math.pow(M2, 2) - Math.pow(A3, 2) * Math.pow(M3, 2) + Math.pow(M1, 2) * Math.pow(a1, 2) + Math.pow(M2, 2) * Math.pow(a2, 2) + Math.pow(M3, 2) * Math.pow(a3, 2)),
        3 * A1 * A2 * Math.pow(M1, 2) * Math.pow(M2, 2) * Math.pow(Q, 2),
        3 * A1 * A3 * Math.pow(M1, 2) * Math.pow(M3, 2) * Math.pow(Q, 2),
        1
    ]
    const row2 = [
        3 * A2 * Math.pow(M1, 2) * Math.pow(M2, 2) * Math.pow(Q, 2),
        3 * A1 * Math.pow(M1, 2) * Math.pow(M2, 2) * Math.pow(Q, 2),
        0,
        0
    ]
    const row3 = [
        3 * A1 * Math.pow(M1, 2) * Math.pow(M2, 2) * Math.pow(Q, 2),
        9 * A2 * Math.pow(M2, 4) * Math.pow(Q, 2),
        3 * A3 * Math.pow(M2, 2) * Math.pow(M3, 2) * Math.pow(Q, 2),
        0
    ]
    const row4 = [
        0,
        3 * A3 * Math.pow(M2, 2) * Math.pow(M3, 2) * Math.pow(Q, 2),
        3 * A2 * Math.pow(M2, 2) * Math.pow(M3, 2) * Math.pow(Q, 2),
        0
    ]
    const row5 = [
        3 * A1 * A3 * Math.pow(M1, 2) * Math.pow(M3, 2) * Math.pow(Q, 2),
        3 * A2 * A3 * Math.pow(M2, 2) * Math.pow(M3, 2) * Math.pow(Q, 2),
        3 * Math.pow(A3, 2) * Math.pow(M3, 4) * Math.pow(Q, 2) + (1 / 2) * Math.pow(M3, 4) - 3 / 2 * Math.pow(M3, 2) * Math.pow(Q, 2) * (-Math.pow(A1, 2) * Math.pow(M1, 2) - Math.pow(A2, 2) * Math.pow(M2, 2) - Math.pow(A3, 2) * Math.pow(M3, 2) + Math.pow(M1, 2) * Math.pow(a1, 2) + Math.pow(M2, 2) * Math.pow(a2, 2) + Math.pow(M3, 2) * Math.pow(a3, 2)),
        -1
    ]
    const jacobian = [
        row1,
        row2,
        row3,
        row4,
        row5
    ]

    return [equation, jacobian]
}

function calculateForces(hCoefficients, wCoefficients, Q, boundary) {
    const [a1, a2, a3] = hCoefficients;
    const [A1, A2, A3] = wCoefficients;

    const SResult = newtonRaphson(FSymmetric, [a1, a2, a3, -50, 50], [Q, a1, a2, a3]);
    const SForce = tempF(a1, a2, a3, SResult[0], SResult[1], SResult[2], Q);
    const STravel = calculateYCoordinate(hCoefficients, boundary, 0.5) - calculateYCoordinate(SResult.slice(0, 3), boundary, 0.5);

    const SBResult = newtonRaphson(FSymmetric, [A1, A2, A3, -50, -50], [Q, a1, a2, a3]);
    const SBForce = tempF(a1, a2, a3, SBResult[0], SBResult[1], SBResult[2], Q);
    const SBTravel = calculateYCoordinate(hCoefficients, boundary, 0.5) - calculateYCoordinate(SBResult.slice(0, 3), boundary, 0.5);

    const ASResult = newtonRaphson(FAsymmetric, [a1, a2, a3, 0], [Q, a1, a2, a3]);
    const ASForce = ASResult[3];
    const ASTravel = calculateYCoordinate(hCoefficients, boundary, 0.5) - calculateYCoordinate(ASResult.slice(0, 3), boundary, 0.5);

    const ASBResult = newtonRaphson(FAsymmetric, [A1, A2, A3, 0], [Q, a1, a2, a3]);
    const ASBForce = ASBResult[3];
    const ASBTravel = calculateYCoordinate(hCoefficients, boundary, 0.5) - calculateYCoordinate(ASBResult.slice(0, 3), boundary, 0.5);

    return { x: [ASTravel, STravel, SBTravel, ASBTravel], y: [ASForce, SForce, SBForce, ASBForce] }
}


function tempF(a1, a2, a3, A1, A2, A3, Q) {
    const M1 = PI;
    const M2 = 2 * PI;
    const M3 = 3 * PI;
    const F = (1 / 2) * pow(M3, 2) * (3 * pow(A1, 2) * A3 * pow(M1, 2) * pow(Q, 2) + 3 * pow(A2, 2) * A3 * pow(M2, 2) * pow(Q, 2) + 3 * pow(A3, 3) * pow(M3, 2) * pow(Q, 2) - 3 * A3 * pow(M1, 2) * pow(Q, 2) * pow(a1, 2) - 3 * A3 * pow(M2, 2) * pow(Q, 2) * pow(a2, 2) - 3 * A3 * pow(M3, 2) * pow(Q, 2) * pow(a3, 2) + A3 * pow(M3, 2) - pow(M3, 2) * a3);
    return F;
}

function newtonRaphson(equationJacobian, x0, parameters) {
    const tolerance = 0.001;
    const maxIteration = 100;
    let x = x0;
    for (let i = 0; i < maxIteration; i++) {
        const array = equationJacobian(x, parameters);
        const equation = array[0];
        const jacobian = array[1];
        if (norm(equation) < tolerance) {
            console.log("Converged!", i + 1);
            break;
        }
        else {
            const inverse = math.pinv(jacobian);
            x = subtract(x, multiply(inverse, equation));
        }
    }
    return x;
}


function allPositive(eigenvalues) {
    return eigenvalues.every((eigenvalue) => eigenvalue >= 0);
}


function bisect(hCoefficients, Q, boundary, a0, b0, maxIteration, tolerance) {
    let a = a0;
    let b = b0;
    for (let i = 0; i < maxIteration; i++) {
        let c = (a + b) / 2;
        const fC = A1Equation(c, hCoefficients, Q, boundary);
        //console.log(abs(a - b), fC);
        if (abs(fC) < tolerance) {
            //console.log(fC);
            console.log("Converged!", c, fC);
            return c;
        }
        else {
            const fA = A1Equation(a, hCoefficients, Q, boundary);
            if (fA * fC < 0) {
                b = c;
            }
            else a = c;
        }
    }
}


function linspace(start, end, number) {
    const result = new Array(number);
    for (let i = 0; i < number; i++) {
        result[i] = ((number - 1) * start + i * (end - start)) / (number - 1);
    }
    return result;
}

function multiply(A, b) {
    const length = A.length;
    const result = new Array(length);
    for (let i = 0; i < length; i++) {
        let sum = 0;
        for (let j = 0; j < b.length; j++) {
            sum += A[i][j] * b[j];
        }
        result[i] = sum;
    }
    return result;
}

function subtract(b, a) {
    const length = b.length;
    const result = new Array(length);
    for (let i = 0; i < length; i++) {
        result[i] = b[i] - a[i];
    }
    return result;
}

function norm(vector) {
    let sum = 0;
    vector.forEach((component) => sum += component * component);
    return Math.sqrt(sum);
}

function inRange(number, start, end) {
    return number >= start && number <= end;
}