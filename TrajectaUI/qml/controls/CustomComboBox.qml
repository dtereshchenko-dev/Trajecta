pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Controls.Basic
import TrajectaUI
import "../shared/UI.js" as UI

ComboBox {
    id: control
    implicitHeight: UI.controlHeight
    implicitWidth: 190
    leftPadding: 10
    rightPadding: 30
    opacity: enabled ? 1 : 0.45
    InteractionCursor {}
    contentItem: CustomText { text: control.displayText; font.pixelSize: 12; verticalAlignment: Text.AlignVCenter }
    indicator: CustomIcon { name: "chevron-right"; rotation: 90; size: 13; x: control.width - 23; y: (control.height - height) / 2; color: AppTheme.muted }
    background: Rectangle { color: AppTheme.base; radius: UI.radius; border.color: control.activeFocus ? AppTheme.accent : control.hovered ? AppTheme.controlLineHover : AppTheme.controlLine }
    delegate: ItemDelegate {
        id: option
        required property int index
        required property var modelData
        width: control.width
        height: UI.controlHeight
        highlighted: control.highlightedIndex === index
        contentItem: CustomText { text: control.textRole ? option.modelData[control.textRole] : option.modelData; font.pixelSize: 12 }
        background: Rectangle { color: option.highlighted ? AppTheme.hover : AppTheme.surface }
        InteractionCursor {}
    }
    popup: Popup {
        y: control.height + 4
        width: control.width
        padding: 1
        implicitHeight: Math.min(contentItem.implicitHeight + 2, 280)
        background: Rectangle { color: AppTheme.surface; border.color: AppTheme.controlLine; radius: UI.radius }
        contentItem: ListView {
            clip: true
            implicitHeight: contentHeight
            model: control.popup.visible ? control.delegateModel : null
            currentIndex: control.highlightedIndex
            ScrollBar.vertical: CustomScrollBar {}
        }
    }
}
