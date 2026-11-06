import QtQuick
import QtQuick.Layouts
import TrajectaUI

RowLayout {
    id: filters
    required property ObservationTableModel observations
    spacing: 12
    CustomTextField {
        objectName: "observations-filter"
        Layout.fillWidth: true
        Layout.maximumWidth: 380
        Layout.minimumWidth: 180
        placeholderText: qsTr("Date, station, program, catalogue")
        Accessible.name: qsTr("Search observations")
        leftPadding: 33
        rightPadding: 31
        text: filters.observations.filterText
        onTextEdited: filters.observations.filterText = text
        CustomIcon { anchors.left: parent.left; anchors.leftMargin: 10; anchors.verticalCenter: parent.verticalCenter; name: "search"; size: 14; color: AppTheme.muted }
        CustomButton {
            objectName: "observations-clear-search"
            anchors.right: parent.right
            anchors.verticalCenter: parent.verticalCenter
            visible: parent.text.length > 0
            iconName: "x"
            flat: true
            width: 30
            height: 30
            accessibleLabel: qsTr("Clear search")
            onClicked: filters.observations.filterText = ""
        }
    }
    CustomComboBox {
        objectName: "observations-station-filter"
        Layout.preferredWidth: 148
        model: [qsTr("All stations")].concat(filters.observations.stationResiduals.stations)
        currentIndex: filters.observations.stationFilter.length ? filters.observations.stationResiduals.stations.indexOf(filters.observations.stationFilter) + 1 : 0
        onActivated: filters.observations.stationFilter = currentIndex > 0 ? filters.observations.stationResiduals.stations[currentIndex - 1] : ""
        Accessible.name: qsTr("Filter observations by station")
    }
    CustomCheckBox {
        objectName: "observations-excluded-filter"
        text: qsTr("Excluded only")
        checked: filters.observations.excludedOnly
        onToggled: filters.observations.excludedOnly = checked
    }
    Item { Layout.fillWidth: true }
    CustomText { text: qsTr("%1 / %2 rows").arg(filters.observations.count).arg(filters.observations.totalCount); secondary: true; font.pixelSize: 11 }
}
