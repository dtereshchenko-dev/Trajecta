import QtQuick
import QtQuick.Layouts
import TrajectaUI

ColumnLayout {
    id: section
    objectName: "station-diagnostics"
    required property StationResidualModel statistics
    signal observationsRequested(string station)
    spacing: 0
    Rectangle { Layout.fillWidth: true; implicitHeight: 1; color: AppTheme.line; Layout.bottomMargin: 24 }
    RowLayout {
        Layout.fillWidth: true
        Layout.bottomMargin: 20
        spacing: 24
        ColumnLayout {
            Layout.fillWidth: false
            spacing: 6
            CustomText { text: qsTr("Residuals by station"); font.pixelSize: 14; font.weight: Font.DemiBold }
            CustomText { text: qsTr("%1 stations \u00b7 %2 included \u00b7 %3 excluded").arg(section.statistics.count).arg(section.statistics.included).arg(section.statistics.excluded); secondary: true; font.pixelSize: 11 }
        }
        Item { Layout.fillWidth: true }
        CustomText { text: qsTr("Sort by"); secondary: true; font.pixelSize: 11 }
        CustomComboBox {
            objectName: "station-diagnostics-sort"
            Layout.preferredWidth: 148
            model: [qsTr("Included count"), qsTr("RMS"), qsTr("Mean offset"), qsTr("Station code")]
            currentIndex: section.statistics.sortMode
            onActivated: section.statistics.sortMode = currentIndex
            Accessible.name: qsTr("Sort station diagnostics")
        }
        CustomTextLink { objectName: "all-station-observations"; Layout.fillWidth: false; text: qsTr("All observations"); iconName: "arrow-up-right"; onClicked: section.observationsRequested("") }
    }
    RowLayout {
        Layout.fillWidth: true
        Layout.bottomMargin: 10
        CustomText { text: qsTr("Included residuals \u00b7 arcsec"); secondary: true; font.pixelSize: 10; Layout.fillWidth: true }
        CustomText { text: qsTr("Unweighted RMS per coordinate"); secondary: true; font.pixelSize: 10 }
    }
    StationResidualTable { Layout.fillWidth: true; statistics: section.statistics; onStationRequested: function(station) { section.observationsRequested(station) } }
    CustomText {
        visible: section.statistics.missing > 0
        Layout.fillWidth: true
        Layout.topMargin: 12
        text: qsTr("%1 included observations without finite residuals; statistics use the available pairs.").arg(section.statistics.missing)
        secondary: true
        font.pixelSize: 11
        wrapMode: Text.Wrap
    }
}
