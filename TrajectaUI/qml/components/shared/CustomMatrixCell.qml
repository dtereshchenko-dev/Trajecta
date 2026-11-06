import QtQuick
import TrajectaUI
import "../../shared/UI.js" as UI

CustomButton {
    id: cell
    property var coefficient
    implicitWidth: 90
    implicitHeight: 68
    padding: 6
    font.family: UI.monoFont
    font.pixelSize: 11
    font.weight: Font.Medium
    contentItem: CustomText {
        objectName: cell.objectName + "-text"
        text: cell.text
        font: cell.font
        horizontalAlignment: Text.AlignHCenter
    }
    background: Rectangle {
        radius: 3
        color: Number.isFinite(cell.coefficient)
            ? UI.mix(cell.coefficient < 0 ? AppTheme.blue : AppTheme.accent, AppTheme.base, 0.05 + Math.abs(cell.coefficient) * 0.33)
            : AppTheme.well
        border.width: cell.checked || cell.hovered || cell.activeFocus ? 2 : 0
        border.color: cell.checked ? AppTheme.ink : AppTheme.accent
    }
}
