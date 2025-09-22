importScripts("Lib/prism.js");

onmessage = function (e) {
    const result = e.data;
    let text = `fabricatedProfile = {\n    \"coefficients\": ${arrayToString(result.fabricatedProfile.coefficients)},\n    \"coordinates\": ${coordinateArrayToString(result.fabricatedProfile.coordinates)}\n}`;
    result.toggledProfiles.forEach((profile, index) => {
        text += "\n"
        let stable = "False";
        if (profile.stable) {
            stable = "True";
        }
        text += `toggledProfile${index + 1} = {\n    \"coefficients\": ${arrayToString(profile.coefficients)},\n    \"coordinates\": ${coordinateArrayToString(profile.coordinates)},\n    \"stable\": ${stable}\n}`;
    })
    postMessage(Prism.highlight(text, Prism.languages.python, "python"));
};

function arrayToString(array) {
    let string = "[";
    array.forEach((number, index) => {
        string += number.toString();
        if (index < array.length - 1) {
            string += ", "
        }
    })
    string += "]";
    return string;
}

function coordinateArrayToString(array) {
    let string = "[";
    array.forEach((pair, index) => {
        string += "(" + pair[0].toString() + ", " + pair[1].toString() + ")";
        if (index < array.length - 1) {
            string += ", "
        }
    })
    string += "]";
    return string;
}