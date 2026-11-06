pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

ColumnLayout {
    id: form
    property string fieldPrefix: "search"
    property string seedUtc: ""
    property string sourceIdentity: ""
    property string lastSourceIdentity: ""
    property bool allowEventWindow: false
    property int intervalMode: allowEventWindow ? 0 : 1
    property string startUtc: ""
    property string endUtc: ""
    property string centerUtc: ""
    property real halfWindowHours: 24
    property bool initialized: false
    readonly property real yearMilliseconds: 365.25 * 86400000
    readonly property real durationYears: (parseUtc(endUtc) - parseUtc(startUtc)) / yearMilliseconds
    readonly property string validationError: {
        if (intervalMode === 0) {
            if (centerUtc.length && !Number.isFinite(parseUtc(centerUtc))) return qsTr("Invalid window centre.")
            if (!Number.isFinite(halfWindowHours) || halfWindowHours <= 0) return qsTr("Half-window must be positive.")
        } else {
            if (startUtc.length && !Number.isFinite(parseUtc(startUtc))) return qsTr("Invalid start date.")
            if (endUtc.length && !Number.isFinite(parseUtc(endUtc))) return qsTr("Invalid end date.")
            if (Number.isFinite(durationYears) && durationYears <= 0) return qsTr("End must be after start.")
        }
        return ""
    }
    spacing: 18

    function parseUtc(text) {
        const parts = /^(\d{4})-(\d{2})-(\d{2})[ T](\d{2}):(\d{2}):(\d{2})(?:\.(\d{1,3}))?$/.exec(text.trim())
        if (!parts || Number(parts[1]) < 1) return NaN
        const date = new Date(0)
        date.setUTCFullYear(Number(parts[1]), Number(parts[2]) - 1, Number(parts[3]))
        date.setUTCHours(Number(parts[4]), Number(parts[5]), Number(parts[6]), Number((parts[7] || "").padEnd(3, "0")))
        if (date.getUTCFullYear() !== Number(parts[1]) || date.getUTCMonth() + 1 !== Number(parts[2])
                || date.getUTCDate() !== Number(parts[3]) || date.getUTCHours() !== Number(parts[4])
                || date.getUTCMinutes() !== Number(parts[5]) || date.getUTCSeconds() !== Number(parts[6])) return NaN
        return date.getTime()
    }
    function formatUtc(milliseconds) {
        const date = new Date(milliseconds)
        if (!Number.isFinite(milliseconds) || !Number.isFinite(date.getTime()) || date.getUTCFullYear() < 1 || date.getUTCFullYear() > 9999) return ""
        return date.toISOString().replace("T", " ").replace("Z", "")
    }
    function initialize() {
        if (initialized || !Number.isFinite(parseUtc(seedUtc))) return
        startUtc = seedUtc
        endUtc = formatUtc(parseUtc(seedUtc) + 30 * yearMilliseconds)
        centerUtc = formatUtc(parseUtc(seedUtc) + halfWindowHours * 3600000)
        initialized = true
    }
    function setYears(years) {
        initialized = true
        endUtc = Number.isFinite(years) && years > 0 ? formatUtc(parseUtc(startUtc) + years * yearMilliseconds) : ""
    }
    onSeedUtcChanged: initialize()
    onSourceIdentityChanged: {
        if (!sourceIdentity.length) return
        if (lastSourceIdentity.length && lastSourceIdentity !== sourceIdentity) {
            initialized = false
            startUtc = ""
            endUtc = ""
            centerUtc = ""
            Qt.callLater(form.initialize)
        }
        lastSourceIdentity = sourceIdentity
    }
    Component.onCompleted: initialize()

    GridLayout {
        Layout.fillWidth: true
        columns: width < 460 ? 1 : 2
        columnSpacing: 20
        rowSpacing: 18
        FormField {
            visible: form.allowEventWindow
            label: qsTr("Event date")
            CustomComboBox {
                objectName: form.fieldPrefix + "-date-mode"
                Layout.fillWidth: true
                model: [qsTr("Set date"), qsTr("Date range")]
                currentIndex: form.intervalMode
                Accessible.name: qsTr("Event date mode")
                onActivated: form.intervalMode = currentIndex
            }
        }
        FormField {
            visible: form.intervalMode === 0
            label: qsTr("Window centre, UTC")
            CustomTextField {
                objectName: form.fieldPrefix + "-center"
                Layout.fillWidth: true
                text: form.centerUtc
                placeholderText: "YYYY-MM-DD HH:MM:SS"
                font.family: UI.monoFont
                font.pixelSize: 12
                Accessible.name: qsTr("Window centre UTC")
                onTextEdited: { form.initialized = true; form.centerUtc = text }
            }
        }
        FormField {
            visible: form.intervalMode === 0
            label: qsTr("Half-window, hours")
            CustomNumberField {
                objectName: form.fieldPrefix + "-half-window"
                Layout.fillWidth: true
                value: form.halfWindowHours
                minimumValue: 0
                exclusiveMinimum: true
                Accessible.name: qsTr("Half-window in hours")
                onValueEdited: value => form.halfWindowHours = value
            }
        }
        FormField {
            visible: form.intervalMode === 1
            label: qsTr("Start, UTC")
            CustomTextField {
                objectName: form.fieldPrefix + "-start"
                Layout.fillWidth: true
                text: form.startUtc
                placeholderText: "YYYY-MM-DD HH:MM:SS"
                font.family: UI.monoFont
                font.pixelSize: 12
                Accessible.name: qsTr("Start UTC")
                onTextEdited: { form.initialized = true; form.startUtc = text }
            }
        }
        FormField {
            visible: form.intervalMode === 1
            label: qsTr("End, UTC")
            CustomTextField {
                objectName: form.fieldPrefix + "-end"
                Layout.fillWidth: true
                text: form.endUtc
                placeholderText: "YYYY-MM-DD HH:MM:SS"
                font.family: UI.monoFont
                font.pixelSize: 12
                Accessible.name: qsTr("End UTC")
                onTextEdited: { form.initialized = true; form.endUtc = text }
            }
        }
        FormField {
            visible: form.intervalMode === 1
            label: qsTr("Duration, years")
            CustomNumberField {
                objectName: form.fieldPrefix + "-years"
                Layout.fillWidth: true
                value: Number.isFinite(form.durationYears) ? Number(form.durationYears.toFixed(6)) : NaN
                minimumValue: 0
                exclusiveMinimum: true
                Accessible.name: qsTr("Duration in Julian years")
                onValueEdited: value => form.setYears(value)
            }
        }
    }
    CustomText {
        Layout.fillWidth: true
        visible: text.length > 0
        text: form.validationError
        color: AppTheme.amber
        font.pixelSize: 11
        wrapMode: Text.Wrap
    }
}
