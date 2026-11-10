/* The browser handles navigation and controls only. Python owns every cleaning
   calculation and draws the original Matplotlib figure inside the iframe. */
"use strict";

const chooseDirectoryButton = document.getElementById("choose-directory");
const selectedDirectory = document.getElementById("selected-directory");
const subjectSelect = document.getElementById("subject");
const activitySelect = document.getElementById("activity");
const recordingSelect = document.getElementById("recording");
const parameterForm = document.getElementById("parameter-form");
const parameterFields = document.getElementById("parameter-fields");
const applyButton = document.getElementById("apply-button");
const resetButton = document.getElementById("reset-button");
const plotFrame = document.getElementById("plot-frame");
const statusText = document.getElementById("status");
const previewState = document.getElementById("preview-state");

let recordingTree = {};
let parameterLabels = {};
let savedParameters = null;
let requestInProgress = false;

function showStatus(message, isError = false) {
    // Text content also keeps filenames and server errors out of HTML markup.
    statusText.textContent = message;
    statusText.hidden = !message;
    statusText.classList.toggle("error", isError);
}

function setBusy(isBusy) {
    // Disable navigation while a request runs so a late response cannot attach
    // one recording's settings or figure to a different selected recording.
    requestInProgress = isBusy;
    chooseDirectoryButton.disabled = isBusy;
    for (const select of [subjectSelect, activitySelect, recordingSelect]) {
        select.disabled = isBusy || select.options.length === 0;
    }
    parameterFields.disabled = isBusy || savedParameters === null;
    applyButton.disabled = isBusy || savedParameters === null;
    resetButton.disabled = isBusy || savedParameters === null;
    applyButton.textContent = isBusy ? "Calculating…" : "Apply & replot";
}

async function postRequest(path, body) {
    const response = await fetch(path, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(body),
    });
    const result = await response.json();
    if (!response.ok) throw new Error(result.error || `Request failed (${response.status}).`);
    return result;
}

function fillSelect(select, values) {
    select.replaceChildren(...values.map(value => new Option(value, value)));
    select.disabled = values.length === 0;
}

function fillActivities() {
    fillSelect(activitySelect, Object.keys(recordingTree[subjectSelect.value] || {}));
    fillRecordings();
}

function fillRecordings() {
    fillSelect(recordingSelect, recordingTree[subjectSelect.value]?.[activitySelect.value] || []);
}

// Group related controls without changing the metadata names sent to Python.
// Bounds retain their endpoints; spike detection uses a strict > comparison.
const parameterGroups = [
    { title: "Pupil diameter (mm)", note: "Discard pupil measurements below Low or above High. Values exactly at either bound are kept.",
      fields: [["pupil_diameter_thresholds", "Low", 0], ["pupil_diameter_thresholds", "High", 1]] },
    { title: "Sudden eye changes (mm/sample)", note: "Mark sudden changes between consecutive samples as spikes when larger than these limits. Equality is allowed.",
      fields: [["pupil_change_threshold", "Pupil"], ["aperture_change_threshold", "Mean aperture"]] },
    { title: "Gaze position bounds (px)", note: "Discard gaze outside the image: below 0 or above the X/Y limit. Both endpoints are kept.",
      fields: [["gaze_x_upper_bound", "X maximum"], ["gaze_y_upper_bound", "Y maximum"]] },
    { title: "Sudden gaze movement (px/sample)", note: "Flag gaze that moves farther than this distance between samples. Equality is allowed; the cleaner combines this flag with blink detection.",
      fields: [["gaze_spike_threshold", "Distance threshold"]] },
    { title: "Blink padding (ns)", note: "Extend each blink by these amounts before and after it. Samples at both extended endpoints are included. 1,000,000 ns = 1 ms.",
      fields: [["blink_pre_ns", "Before"], ["blink_post_ns", "After"]] },
    { title: "Fill limits (sample-index distance)", note: "Fill missing stretches only when their sample-index distance is at most this limit. The distance includes the bounding valid samples.",
      fields: [["general_max_gap", "Other gaps"], ["blink_max_gap", "Blink gaps"]] },
    { title: "Window sizes (samples)", note: "Mask width extends rejected regions to neighboring samples. Smoothing width reduces short fluctuations. Both are at least 1; smoothing must be odd.",
      fields: [["dilation_samples", "Mask width"], ["smoothing_samples", "Smoothing width"]] },
];

function displayParameters(parameters) {
    // Every input has a visible label, including separate low/high bounds.
    // The data-parameter attribute preserves the existing request format.
    parameterFields.replaceChildren();
    const integerParameters = ["blink_pre_ns", "blink_post_ns", "gaze_x_upper_bound",
        "gaze_y_upper_bound", "dilation_samples", "general_max_gap", "blink_max_gap"];
    for (const group of parameterGroups) {
        const section = document.createElement("div");
        section.className = "parameter-group";
        const title = document.createElement("h3");
        title.textContent = group.title;
        section.append(title);
        const inputs = document.createElement("div");
        inputs.className = group.fields.length === 2 ? "parameter-pair" : "parameter-single";
        for (const [name, labelText, index = 0] of group.fields) {
            const field = document.createElement("div");
            const label = document.createElement("label");
            label.textContent = labelText;
            label.htmlFor = `${name}-${index}`;
            const input = document.createElement("input");
            input.id = label.htmlFor;
            input.type = "number";
            input.required = true;
            input.min = ["dilation_samples", "smoothing_samples"].includes(name) ? "1" : "0";
            input.step = name === "smoothing_samples" ? "2" :
                (integerParameters.includes(name) ? "1" : "any");
            input.value = Array.isArray(parameters[name]) ? parameters[name][index] : parameters[name];
            input.dataset.parameter = name;
            input.setAttribute("aria-label", `${group.title}: ${labelText}`);
            field.append(label, input);
            inputs.append(field);
        }
        const note = document.createElement("p");
        note.className = "parameter-note";
        note.textContent = group.note;
        section.append(inputs, note);
        parameterFields.append(section);
    }
}

function collectParameters() {
    // Preserve the pupil-bound pair while sending all other fields as scalars.
    const parameters = {};
    for (const name of Object.keys(parameterLabels)) {
        const values = [...parameterFields.querySelectorAll(`[data-parameter="${name}"]`)]
            .map(input => Number(input.value));
        parameters[name] = values.length === 2 ? values : values[0];
    }
    return parameters;
}

async function showRecording(parameters = null) {
    if (requestInProgress || !recordingSelect.value) return;
    const selection = {
        subject: subjectSelect.value,
        activity: activitySelect.value,
        recording: recordingSelect.value,
        parameters,
    };
    // A new recording must never show the previous recording's plot or controls.
    if (parameters === null) {
        savedParameters = null;
        parameterFields.replaceChildren();
        plotFrame.hidden = true;
        plotFrame.removeAttribute("src");
        document.getElementById("empty-state").hidden = false;
        previewState.textContent = "Loading recording";
    }
    previewState.hidden = false;
    setBusy(true);
    document.getElementById("recording-title").textContent =
        `${selection.subject} / ${selection.activity} / ${selection.recording}`;
    showStatus("Running the Python cleaner and preparing the comparison…");
    try {
        const result = await postRequest("/api/plot", selection);
        savedParameters = result.saved_parameters;
        displayParameters(result.parameters);
        plotFrame.src = result.figure_url;
        plotFrame.hidden = false;
        document.getElementById("empty-state").hidden = true;
        const changed = Object.keys(savedParameters).some(name =>
            JSON.stringify(savedParameters[name]) !== JSON.stringify(result.parameters[name]));
        previewState.textContent = changed ? "Adjusted preview" : "Saved settings";
        showStatus("");
    } catch (error) {
        if (savedParameters === null) previewState.textContent = "Unable to load";
        showStatus(error.message, true);
    } finally {
        setBusy(false);
    }
}

chooseDirectoryButton.addEventListener("click", async () => {
    if (requestInProgress) return;
    setBusy(true);
    showStatus("Choose a folder in the system dialog…");
    let directoryLoaded = false;
    try {
        const result = await postRequest("/api/choose-directory", {});
        if (result.cancelled) {
            showStatus("");
            return;
        }
        selectedDirectory.textContent = result.directory;
        selectedDirectory.hidden = false;
        recordingTree = result.subjects;
        fillSelect(subjectSelect, Object.keys(recordingTree));
        fillActivities();
        directoryLoaded = true;
    } catch (error) {
        showStatus(error.message, true);
    } finally {
        setBusy(false);
    }
    if (directoryLoaded) await showRecording();
});

subjectSelect.addEventListener("change", () => { fillActivities(); showRecording(); });
activitySelect.addEventListener("change", () => { fillRecordings(); showRecording(); });
recordingSelect.addEventListener("change", () => showRecording());
parameterForm.addEventListener("submit", event => {
    event.preventDefault();
    if (parameterForm.reportValidity()) showRecording(collectParameters());
});
resetButton.addEventListener("click", () => showRecording(savedParameters));
parameterForm.addEventListener("input", () => {
    previewState.textContent = "Unapplied changes";
});

// Read optional startup configuration, then let the user choose when to load.
fetch("/api/config").then(response => {
    if (!response.ok) throw new Error("Unable to read viewer configuration.");
    return response.json();
}).then(configuration => {
    parameterLabels = configuration.labels;
}).catch(error => showStatus(error.message, true));
