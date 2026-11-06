pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

ColumnLayout {
    id: details
    objectName: "covariance-details"
    required property CovarianceTableModel matrix
    readonly property var pair: matrix.selection
    spacing: 0
    function scientific(value, digits) { return Number.isFinite(value) ? value.toExponential(digits) : "\u2014" }
    function sigmaLabel(label, unit) { return "\u03c3 " + (label || "") + (unit ? " [" + unit + "]" : "") }
    CustomText {
        objectName: "covariance-pair-title"
        text: details.pair.rowLabel ? details.pair.rowLabel + " / " + details.pair.columnLabel : "\u2014"
        font.pixelSize: 14
        font.weight: Font.DemiBold
        Layout.fillWidth: true
        Layout.bottomMargin: 20
    }
    KeyValueList {
        Layout.fillWidth: true
        rows: [
            [qsTr("Correlation"), UI.number(details.pair.correlation, 8)],
            [qsTr("Covariance"), details.scientific(details.pair.covariance, 7)],
            [details.sigmaLabel(details.pair.rowLabel, details.pair.rowUnit), details.scientific(details.pair.rowSigma, 6)],
            [details.sigmaLabel(details.pair.columnLabel, details.pair.columnUnit), details.scientific(details.pair.columnSigma, 6)]
        ]
    }
    Rectangle { Layout.fillWidth: true; implicitHeight: 1; color: AppTheme.line; Layout.topMargin: 26; Layout.bottomMargin: 24 }
    CustomText { text: qsTr("Weak direction"); font.pixelSize: 14; font.weight: Font.DemiBold; Layout.bottomMargin: 20 }
    KeyValueList {
        Layout.fillWidth: true
        visible: details.matrix.weakDirection.length > 0
        rows: [[qsTr("Weak \u03c3 (EIG)"), details.scientific(details.matrix.weakSigma, 5)]].concat(
            details.matrix.weakDirection.map(function(row) { return [row[0], UI.signedNumber(row[1], 5)] }))
    }
    CustomText {
        visible: details.matrix.weakDirection.length === 0
        text: qsTr("Weak direction unavailable")
        secondary: true
        font.pixelSize: 11
        Layout.fillWidth: true
        wrapMode: Text.Wrap
    }
    CustomText {
        text: details.matrix.equinoctial ? qsTr("EIG / WEA \u00b7 EQU, \u03bb in radians") : qsTr("EIG / WEA \u00b7 CAR, AU and AU/day")
        secondary: true
        font.pixelSize: 10
        Layout.topMargin: 12
        Layout.fillWidth: true
        wrapMode: Text.Wrap
    }
}
