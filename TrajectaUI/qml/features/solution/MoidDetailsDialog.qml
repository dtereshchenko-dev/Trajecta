pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

BaseDialog {
    id: dialog
    objectName: "moid-details-dialog"
    required property SolutionViewModel solution
    title: qsTr("Minimum orbit distances")
    width: Math.min(550, parent.width - 40)
    contentItem: ColumnLayout {
        spacing: 14
        CustomText { text: dialog.solution.epochMjd + " MJD TDB"; secondary: true; font.pixelSize: 11; Layout.fillWidth: true }
        RowLayout {
            Layout.fillWidth: true
            CustomText { text: qsTr("Planet"); secondary: true; Layout.fillWidth: true }
            CustomText { text: "AU"; secondary: true; Layout.preferredWidth: 140; horizontalAlignment: Text.AlignRight }
            CustomText { text: "km"; secondary: true; Layout.preferredWidth: 130; horizontalAlignment: Text.AlignRight }
        }
        Rectangle { Layout.fillWidth: true; implicitHeight: 1; color: AppTheme.line }
        Repeater {
            model: dialog.solution.planetMoids
            RowLayout {
                id: row
                required property var modelData
                Layout.fillWidth: true
                CustomText { text: row.modelData.planet; Layout.fillWidth: true; font.pixelSize: 12 }
                CustomText { text: UI.number(row.modelData.au, 10); numeric: true; font.pixelSize: 11; Layout.preferredWidth: 140; horizontalAlignment: Text.AlignRight }
                CustomText { text: UI.number(row.modelData.km, 1); numeric: true; font.pixelSize: 11; Layout.preferredWidth: 130; horizontalAlignment: Text.AlignRight }
            }
        }
    }
}
