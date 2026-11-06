import QtQuick
import TrajectaUI

Column {
    id: section
    required property string title
    required property var rows
    spacing: 14
    Rectangle { width: parent.width; height: 1; color: AppTheme.line }
    CustomText { text: section.title; font.pixelSize: 12; font.weight: Font.DemiBold }
    KeyValueList { width: parent.width; rows: section.rows; spacing: 12 }
}
