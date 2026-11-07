pragma Singleton
import QtQuick
import "UI.js" as UI

QtObject {
    readonly property bool windows: Qt.platform.os === "windows"
    property FontLoader textFont: FontLoader { source: windows ? "../../resources/fonts/Inter.ttc" : "../../resources/fonts/SourceSans3-Regular.ttf" }
    property FontLoader textItalicFont: FontLoader { source: windows ? "" : "../../resources/fonts/SourceSans3-It.ttf" }
    property FontLoader textSemiBoldFont: FontLoader { source: windows ? "" : "../../resources/fonts/SourceSans3-Semibold.ttf" }
    property FontLoader textBoldFont: FontLoader { source: windows ? "" : "../../resources/fonts/SourceSans3-Bold.ttf" }
    property FontLoader monoFont: FontLoader { source: windows ? "../../resources/fonts/JetBrainsMono-Regular.ttf" : "../../resources/fonts/CascadiaMono.ttf" }
    property FontLoader monoBoldFont: FontLoader { source: windows ? "../../resources/fonts/JetBrainsMono-Bold.ttf" : "" }
    property string name: "graphite"
    readonly property var colors: UI.themes[name] || UI.themes.graphite
    readonly property color base: colors.base
    readonly property color surface: colors.surface
    readonly property color nav: colors.nav
    readonly property color line: colors.line
    readonly property color ink: colors.ink
    readonly property color muted: colors.muted
    readonly property color accent: colors.accent
    readonly property color accentText: colors["onAccent"]
    readonly property color blue: colors.blue
    readonly property color red: colors.red
    readonly property color amber: colors.amber
    readonly property color surfaceSecondary: UI.mix(surface, base, 0.86)
    readonly property color well: UI.mix(nav, base, 0.72)
    readonly property color hover: UI.mix(nav, ink, 0.88)
    readonly property color rowHover: UI.mix(surface, accent, 0.91)
    readonly property color lineSoft: UI.mix(line, surface, 0.55)
    readonly property color controlLine: UI.mix(line, ink, 0.82)
    readonly property color controlLineHover: UI.mix(line, ink, 0.55)
    readonly property color accentSoft: UI.mix(accent, surface, 0.13)
    readonly property color accentTint: UI.mix(accent, nav, 0.20)
    readonly property color accentStrong: UI.mix(accent, ink, 0.82)
    readonly property color chartGrid: UI.mix(line, surface, 0.72)
    readonly property color chartTick: UI.mix(muted, surface, 0.88)
    readonly property color chartZero: UI.mix(line, muted, 0.45)
    readonly property color shadow: name === "light" ? "#23312f" : name === "dim" ? "#1b2523" : "#000000"
    readonly property color overlay: Qt.alpha(shadow, 0.28)
}
