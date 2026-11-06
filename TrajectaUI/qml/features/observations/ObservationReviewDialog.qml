import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Layouts
import TrajectaUI

BaseDialog {
    id: dialog
    required property ApplicationSession session
    readonly property var review: session.observationReview
    title: qsTr("Review observations")
    width: Math.min(920, parent.width - 40)
    height: Math.min(680, parent.height - 40)
    acceptText: qsTr("Solve")
    cancelText: qsTr("Cancel")
    acceptPrimary: true
    canAccept: review.remainingOpticalCount >= 3
    onOpened: filter.currentIndex = 0
    onAccepted: session.acceptObservationReview()
    onRejected: session.cancelObservationReview()

    contentItem: ColumnLayout {
        spacing: 12
        RowLayout {
            Layout.fillWidth: true
            CustomText { text: qsTr("Issue"); secondary: true }
            CustomComboBox {
                id: filter
                Layout.preferredWidth: 280
                model: [qsTr("All issues"), qsTr("Same time and position"), qsTr("Nearby observations, conflicting positions"),
                    qsTr("Below horizon"), qsTr("Daylight")]
                onActivated: function(index) {
                    dialog.review.filter = index === 1 ? dialog.review.positionFilter
                        : index === 2 ? dialog.review.differentPositionFilter
                        : index === 3 ? dialog.review.belowHorizonFilter
                        : index === 4 ? dialog.review.daylightFilter : 0
                }
            }
            Item { Layout.fillWidth: true }
            CustomText {
                text: qsTr("Shown: %1 / %2   Excluded: %3").arg(dialog.review.shownCount)
                    .arg(dialog.review.totalCount).arg(dialog.review.excludedCount)
                secondary: true
            }
        }
        RowLayout {
            Layout.fillWidth: true
            spacing: 8
            CustomButton {
                text: qsTr("Select all shown")
                onClicked: dialog.review.setAllShownIncluded(true)
            }
            CustomButton {
                text: qsTr("Deselect all shown")
                onClicked: dialog.review.setAllShownIncluded(false)
            }
            Item { Layout.fillWidth: true }
        }
        Rectangle {
            Layout.fillWidth: true
            Layout.fillHeight: true
            color: AppTheme.base
            border.color: AppTheme.line
            clip: true
            ListView {
                id: list
                anchors.fill: parent
                anchors.margins: 1
                model: dialog.review
                clip: true
                boundsBehavior: Flickable.StopAtBounds
                ScrollBar.vertical: CustomScrollBar {
                    parent: list
                    x: list.width - width
                    y: 0
                    height: list.height
                    policy: ScrollBar.AlwaysOn
                    active: true
                    contentItem: Rectangle {
                        implicitWidth: 6
                        radius: 3
                        color: AppTheme.muted
                    }
                }
                delegate: Rectangle {
                    required property int index
                    required property int groupNumber
                    required property string station
                    required property string observationTime
                    required property string coordinates
                    required property string reason
                    required property bool included
                    required property bool canChange
                    width: ListView.view.width
                    height: 62
                    color: groupNumber % 2 ? AppTheme.surface : AppTheme.base
                    RowLayout {
                        anchors.fill: parent
                        anchors.leftMargin: 12
                        anchors.rightMargin: 18
                        spacing: 12
                        CustomCheckBox {
                            checked: included
                            enabled: canChange
                            onClicked: dialog.review.setIncluded(index, checked)
                        }
                        CustomText { text: groupNumber ? "#" + groupNumber : ""; secondary: true; Layout.preferredWidth: 38 }
                        ColumnLayout {
                            Layout.fillWidth: true
                            spacing: 3
                            CustomText {
                                text: station + "   " + observationTime + "   " + coordinates
                                font.pixelSize: 12
                                elide: Text.ElideRight
                                Layout.fillWidth: true
                            }
                            CustomText {
                                text: reason
                                color: AppTheme.amber
                                font.pixelSize: 10
                                elide: Text.ElideRight
                                Layout.fillWidth: true
                            }
                        }
                    }
                    Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 1; color: AppTheme.line }
                }
            }
        }
    }
}
