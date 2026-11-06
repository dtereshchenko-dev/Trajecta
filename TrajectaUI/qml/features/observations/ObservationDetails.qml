import QtQuick
import TrajectaUI

Column {
    id: details
    required property ObservationTableModel observations
    readonly property var selection: observations.selection
    visible: !!selection.number
    spacing: 22
    CustomText {
        width: parent.width
        text: details.selection.included ? qsTr("Included in fit") : details.selection.replaced ? qsTr("Replaced by MPC") : qsTr("Excluded from fit")
        color: details.selection.included ? AppTheme.accent : AppTheme.red
        font.pixelSize: 11
    }
    KeyValueList {
        width: parent.width
        spacing: 12
        rows: details.selection.radar ? [
            [qsTr("Date") + " " + details.selection.timeScale, details.selection.epoch],
            [qsTr("Transmitter"), details.selection.transmitter],
            [qsTr("Receiver"), details.selection.station],
            [qsTr("Frequency, MHz"), details.selection.frequency],
            [qsTr("Return"), details.selection.returnType]
        ] : details.selection.number ? [
            [qsTr("Date") + " " + details.selection.timeScale, details.selection.epoch],
            [qsTr("Station"), details.selection.station],
            [qsTr("RA"), details.selection.ra],
            [qsTr("Dec"), details.selection.dec]
        ] : []
    }
    ObservationDetailSection {
        width: parent.width
        visible: !details.selection.radar
        title: qsTr("Measurement")
        rows: details.selection.number ? [
            [qsTr("Technique"), details.selection.technique],
            [qsTr("Catalogue"), details.selection.catalog],
            [qsTr("MPC program"), details.selection.program || qsTr("Not supplied")]
        ] : []
    }
    ObservationDetailSection {
        width: parent.width
        visible: !details.selection.radar
        title: qsTr("Weights and debias")
        rows: details.selection.number ? [
            [qsTr("Weighting"), details.selection.weighting],
            [qsTr("\u03c3 RA cos \u03b4, \u2033"), details.selection.sigmaRa],
            [qsTr("\u03c3 Dec, \u2033"), details.selection.sigmaDec],
            [qsTr("Applied debias RA, \u2033"), details.selection.debiasEnabled ? details.selection.biasRa : qsTr("Off")],
            [qsTr("Applied debias Dec, \u2033"), details.selection.debiasEnabled ? details.selection.biasDec : qsTr("Off")]
        ] : []
    }
    ObservationDetailSection {
        width: parent.width
        visible: !details.selection.radar
        title: qsTr("Residuals")
        rows: details.selection.number ? [
            [qsTr("\u0394\u03b1 cos \u03b4, \u2033"), details.selection.raResidual],
            [qsTr("\u0394\u03b4, \u2033"), details.selection.decResidual],
            [qsTr("\u03c7 (measurement)"), details.selection.chi],
            [qsTr("\u03c7 (post-fit)"), details.selection.postFitChi]
        ] : []
    }
    Repeater {
        model: details.selection.radarChannels || []
        ObservationDetailSection {
            required property var modelData
            width: details.width
            title: modelData.label
            rows: [[qsTr("Observed"), modelData.observed], [qsTr("Computed"), modelData.computed],
                [qsTr("O-C"), modelData.residual], [qsTr("Sigma"), modelData.sigma]]
        }
    }
    ObservationDetailSection {
        width: parent.width
        visible: !!details.selection.radar
        title: qsTr("Fit statistics")
        rows: details.selection.radar ? [[qsTr("\u03c7 (measurement)"), details.selection.chi], [qsTr("\u03c7 (post-fit)"), details.selection.postFitChi]] : []
    }
}
