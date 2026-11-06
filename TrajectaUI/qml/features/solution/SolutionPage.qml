import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

BasePage {
    id: page
    objectName: "solution-page"
    required property SolutionViewModel solution
    property int currentTab: 0
    property real viewportWidth: 1440
    property bool canRecalculate: false
    signal observationsRequested(string station)
    signal gaiaRequested()
    signal importRequested()
    title: qsTr("Solution")
    padding: 0
    contentItem: ColumnLayout {
        spacing: 0
        Item {
            Layout.fillWidth: true
            implicitHeight: 47
            RowLayout {
                anchors.fill: parent
                anchors.leftMargin: page.viewportWidth <= 1200 ? 20 : 30
                anchors.rightMargin: anchors.leftMargin
                spacing: 24
                CustomTab { text: qsTr("Overview"); checked: page.currentTab === 0; onClicked: page.currentTab = 0 }
                CustomTab { text: qsTr("Covariance"); checked: page.currentTab === 1; enabled: page.solution.available; onClicked: page.currentTab = 1 }
                CustomTab { objectName: "solution-charts-tab"; text: qsTr("Residuals"); checked: page.currentTab === 2; enabled: page.solution.available; onClicked: page.currentTab = 2 }
                Item { Layout.fillWidth: true }
                CustomButton {
                    visible: page.solution.available && !page.solution.hasGaiaFprObservations
                    text: qsTr("Add Gaia FPR and re-solve")
                    iconName: "download"
                    enabled: page.canRecalculate
                    onClicked: page.gaiaRequested()
                }
                CustomText {
                    visible: page.solution.available && page.solution.hasGaiaObservations
                    text: page.solution.gaiaSummary
                    secondary: true
                    font.pixelSize: 11
                }
                CustomText { text: qsTr("Nominal solution"); secondary: true; font.pixelSize: 11 }
            }
            Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 1; color: AppTheme.line }
        }
        StackLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            currentIndex: page.solution.available ? page.currentTab : 3
            Flickable {
                id: overview
                objectName: "solution-overview"
                clip: true
                contentWidth: width
                contentHeight: summary.implicitHeight + 56
                boundsBehavior: Flickable.StopAtBounds
                ScrollBar.vertical: CustomScrollBar {}
                SolutionSummary {
                    id: summary
                    x: UI.pageMargin(page.viewportWidth)
                    y: 24
                    width: overview.width - 2 * x
                    solution: page.solution
                    viewportWidth: page.viewportWidth
                    onObservationsRequested: function(station) { page.observationsRequested(station) }
                }
            }
            Item {
                CovariancePanel {
                    anchors.fill: parent
                    anchors.leftMargin: UI.pageMargin(page.viewportWidth)
                    anchors.rightMargin: anchors.leftMargin
                    anchors.topMargin: 24
                    matrix: page.solution.covariance
                }
            }
            SolutionCharts {
                diagnostics: page.solution.diagnostics
                observations: page.solution.observations
                pageMargin: UI.pageMargin(page.viewportWidth)
            }
            Item {
                ColumnLayout {
                    anchors.centerIn: parent
                    spacing: 16
                    CustomIcon { name: "orbit"; size: 40; color: AppTheme.accent; Layout.alignment: Qt.AlignHCenter }
                    CustomText { text: qsTr("No calculated solution"); secondary: true }
                    CustomButton { text: qsTr("Open observations"); iconName: "folder-open"; Layout.alignment: Qt.AlignHCenter; onClicked: page.importRequested() }
                }
            }
        }
    }
}
