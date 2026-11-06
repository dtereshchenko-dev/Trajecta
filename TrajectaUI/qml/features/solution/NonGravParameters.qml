pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI

ColumnLayout {
    id: root
    required property SolutionViewModel solution
    property bool detailed: false
    readonly property var rows: detailed ? solution.nonGravParameters : solution.nonGravParameters.filter(row => row.used)
    readonly property int parameterWidth: 76
    readonly property int sigmaWidth: 72
    spacing: 0

    function formatValue(value, label, estimated) {
        if (!Number.isFinite(value) || (!estimated && value === 0)) return "\u2014"
        return value === 0 ? "0" : label === "DT" ? value.toFixed(6) : value.toExponential(6)
    }

    RowLayout {
        Layout.fillWidth: true
        Layout.bottomMargin: root.solution.nonGravEnabled ? 10 : 0
        spacing: 12
        SectionHeading { title: qsTr("Non-gravitational"); Layout.fillWidth: !root.detailed }
        CustomText {
            objectName: "non-grav-model"
            text: root.solution.nonGravModel || "\u2014"
            color: root.detailed && root.solution.nonGravEnabled ? AppTheme.accent : AppTheme.muted
            font.pixelSize: 11
            Layout.maximumWidth: root.width / 2
            wrapMode: Text.Wrap
        }
        Item { visible: root.detailed; Layout.fillWidth: true }
    }
    ParameterTable {
        visible: root.detailed && root.solution.nonGravEnabled
        Layout.fillWidth: true
        rowPrefix: "non-grav"
        showSigma: true
        showStatus: true
        showSignificance: true
        columnSpacing: width >= 700 ? 24 : 12
        parameterShare: 0.13
        sigmaShare: 0.17
        significanceShare: 0.15
        unitShare: 0.15
        statusShare: 0.13
        rows: root.detailed && root.solution.nonGravEnabled ? root.rows.map(row => ({
            label: row.label, value: root.formatValue(row.value, row.label, row.estimated), unit: row.unit,
            sigma: Number.isFinite(row.sigma) ? row.sigma.toExponential(3) : "\u2014",
            significance: row.estimated && Number.isFinite(row.value) && Number.isFinite(row.sigma) && row.sigma > 0
                ? (Math.abs(row.value) / row.sigma).toPrecision(3) : "\u2014",
            status: row.status, estimated: row.estimated, used: row.used
        })) : []
    }
    ColumnLayout {
        visible: !root.detailed && root.solution.nonGravEnabled
        Layout.fillWidth: true
        spacing: 0
        RowLayout {
            Layout.fillWidth: true
            Layout.preferredHeight: 28
            spacing: 10
            CustomText { text: qsTr("Parameter"); secondary: true; font.pixelSize: 10; Layout.preferredWidth: root.parameterWidth }
            CustomText { text: qsTr("Value"); secondary: true; font.pixelSize: 10; Layout.fillWidth: true; horizontalAlignment: Text.AlignRight }
            CustomText { text: "\u03c3"; secondary: true; font.pixelSize: 10; Layout.preferredWidth: root.sigmaWidth; horizontalAlignment: Text.AlignRight }
        }
        Repeater {
            model: root.detailed ? [] : root.rows
            Item {
                id: row
                required property var modelData
                objectName: "non-grav-row-" + modelData.label
                Layout.fillWidth: true
                implicitHeight: 34
                Rectangle { width: parent.width; height: 1; color: AppTheme.lineSoft }
                RowLayout {
                    anchors.fill: parent
                    spacing: 10
                    RowLayout {
                        Layout.preferredWidth: root.parameterWidth
                        spacing: 8
                        CustomText { text: row.modelData.label; font.pixelSize: 11; font.weight: Font.DemiBold }
                        CustomText { text: row.modelData.unit; secondary: true; font.pixelSize: 10; Layout.fillWidth: true }
                    }
                    CustomText {
                        objectName: "non-grav-value"
                        text: root.formatValue(row.modelData.value, row.modelData.label, row.modelData.estimated)
                        numeric: true
                        font.pixelSize: 11
                        font.weight: Font.DemiBold
                        Layout.fillWidth: true
                        horizontalAlignment: Text.AlignRight
                    }
                    CustomText {
                        objectName: "non-grav-sigma"
                        text: !row.modelData.estimated ? row.modelData.status
                            : Number.isFinite(row.modelData.sigma) ? row.modelData.sigma.toExponential(3) : "\u2014"
                        numeric: row.modelData.estimated
                        secondary: true
                        font.pixelSize: 11
                        Layout.preferredWidth: root.sigmaWidth
                        horizontalAlignment: Text.AlignRight
                    }
                }
            }
        }
    }
}
