.pragma library

var systemFonts = Qt.fontFamilies()
var font = "Source Sans 3"
var monoFont = "Cascadia Mono"
if (Qt.platform.os === "windows") {
    font = systemFonts.indexOf("Segoe UI Variable Text") >= 0 ? "Segoe UI Variable Text" : systemFonts.indexOf("Segoe UI") >= 0 ? "Segoe UI" : "Inter"
    monoFont = systemFonts.indexOf("Cascadia Mono") >= 0 ? "Cascadia Mono" : systemFonts.indexOf("Consolas") >= 0 ? "Consolas" : "JetBrains Mono"
}
var controlHeight = 33
var rowHeight = 36
var radius = 4
var navigationWidth = 202
var navigationRailWidth = 59
var commandBarHeight = 58
var headingSize = 28
var sectionSize = 14
var bodySize = 13
var captionSize = 11
var themes = Object.freeze({
    graphite: Object.freeze({base: "#16191a", surface: "#1d2122", nav: "#111415", line: "#2e3436", ink: "#dfe5e4", muted: "#8d9998", accent: "#4fae94", onAccent: "#06120f", blue: "#6f9fdc", red: "#e0736a", amber: "#d6a24a"}),
    observatory: Object.freeze({base: "#10151c", surface: "#171d26", nav: "#0b1016", line: "#27303c", ink: "#dbe3ec", muted: "#8794a4", accent: "#46a6b8", onAccent: "#03161a", blue: "#6d9fe0", red: "#e07a74", amber: "#d5a154"}),
    slate: Object.freeze({base: "#262d31", surface: "#2f373b", nav: "#202629", line: "#414a4f", ink: "#dde4e6", muted: "#93a0a4", accent: "#5aae9a", onAccent: "#061512", blue: "#7aa5d8", red: "#dd7c73", amber: "#cfa055"}),
    nocturne: Object.freeze({base: "#000000", surface: "#0c0e0f", nav: "#050607", line: "#232829", ink: "#e6ecea", muted: "#8e9a98", accent: "#57bfa0", onAccent: "#001410", blue: "#74a6e4", red: "#e97b71", amber: "#dba64f"}),
    redlight: Object.freeze({base: "#150e0c", surface: "#1c1412", nav: "#100a08", line: "#332420", ink: "#ecd9cf", muted: "#a48a80", accent: "#cf6a3f", onAccent: "#1a0b04", blue: "#c39a5e", red: "#e2574a", amber: "#d99a3c"}),
    dim: Object.freeze({base: "#e3e8e6", surface: "#eef2f0", nav: "#d9e0dd", line: "#c2ccc8", ink: "#1d2a28", muted: "#586563", accent: "#276558", onAccent: "#ffffff", blue: "#3f6aa5", red: "#a04a43", amber: "#84611f"}),
    light: Object.freeze({base: "#f5f7f7", surface: "#ffffff", nav: "#eef2f1", line: "#dce3e2", ink: "#23312f", muted: "#64736f", accent: "#286e60", onAccent: "#ffffff", blue: "#4877b9", red: "#ae534b", amber: "#9b7026"})
})

function getSVG(path) {
    var icons = ["orbit", "scan-line", "table-2", "folder-open", "folder", "file-text", "play", "square", "chevron-right", "circle-dot", "arrow-up-right", "x", "moon", "sun", "sliders-horizontal", "activity", "cloud-sun", "sigma", "info", "check", "panel-left-close", "panel-left-open", "calendar-range", "download"]
    if (icons.indexOf(path) >= 0 || ["plus", "minus", "rotate-ccw", "view", "pause", "search", "triangle-alert"].indexOf(path) >= 0)
        return "qrc:/qt/qml/TrajectaUI/resources/icons/workspace/" + path + ".svg"
    return ""
}

function pageMargin(width) { return width <= 1200 ? 20 : width >= 1700 ? 40 : 30 }

function mix(a, b, fraction) {
    return Qt.rgba(a.r * fraction + b.r * (1 - fraction), a.g * fraction + b.g * (1 - fraction), a.b * fraction + b.b * (1 - fraction), 1)
}

function number(value, digits) {
    return Number.isFinite(value) ? value.toFixed(digits === undefined ? 6 : digits) : "—"
}

function sigma(value) {
    return Number.isFinite(value) && value > 0 ? value.toExponential(3) : "—"
}

function signedNumber(value, digits) {
    if (!Number.isFinite(value)) return "\u2014"
    var rounded = Number(value.toFixed(digits))
    return (rounded > 0 ? "+" : "") + rounded.toFixed(digits)
}

function groupedNumber(value, digits) {
    if (!Number.isFinite(value)) return "—"
    var parts = number(value, digits).split(".")
    parts[0] = parts[0].replace(/\B(?=(\d{3})+(?!\d))/g, " ")
    return parts.join(".")
}
