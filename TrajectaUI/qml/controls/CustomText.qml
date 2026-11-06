import QtQuick
import TrajectaUI
import "../shared/UI.js" as UI

Text {
    property bool secondary: false
    property bool numeric: false
    color: secondary ? AppTheme.muted : AppTheme.ink
    font.family: numeric ? UI.monoFont : UI.font
    font.pixelSize: UI.bodySize
    font.letterSpacing: 0
    elide: Text.ElideRight
    verticalAlignment: Text.AlignVCenter
    textFormat: Text.PlainText
}
