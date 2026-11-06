pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI

RowLayout {
    id: section
    property string title
    Layout.fillWidth: true
    Layout.topMargin: 10
    spacing: 12
    CustomText { text: section.title; font.pixelSize: 11; font.weight: Font.DemiBold }
    Rectangle { Layout.fillWidth: true; implicitHeight: 1; color: AppTheme.line }
}
