import QtQuick
import TrajectaUI
import "../shared/UI.js" as UI

CustomTextField {
    id: control
    property real value: 0
    property real minimumValue: -Infinity
    property bool exclusiveMinimum: false
    readonly property bool valid: Number.isFinite(value) && (exclusiveMinimum ? value > minimumValue : value >= minimumValue)
    signal valueEdited(real number)
    implicitWidth: 140
    font.family: UI.monoFont
    font.pixelSize: 12
    color: valid ? AppTheme.ink : AppTheme.red
    opacity: enabled ? 1 : 0.45
    function refresh() { text = Number.isFinite(value) ? String(value) : "" }
    onValueChanged: if (!activeFocus) refresh()
    onActiveFocusChanged: if (!activeFocus) refresh()
    Component.onCompleted: refresh()
    onTextEdited: valueEdited(text.trim().length && Number.isFinite(Number(text)) ? Number(text) : NaN)
}
