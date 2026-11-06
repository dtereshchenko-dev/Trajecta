pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import QtQuick3D
import QtQuick3D.Helpers
import TrajectaUI
import "../../shared/UI.js" as UI

Rectangle {
    id: scene
    objectName: "orbit-scene"
    required property SolutionViewModel solution
    property bool showPlanets: true
    property bool showGrid: true
    property bool showLabels: true
    property bool showPerihelion: true
    readonly property OrbitSceneModel sceneModel: solution.scene
    readonly property real gridStep: {
        const desired = Math.max(1e-6, camera.z * 2 * Math.tan(camera.fieldOfView * Math.PI / 360) / 12)
        const magnitude = Math.pow(10, Math.floor(Math.log10(desired)))
        const fraction = desired / magnitude
        return magnitude * (fraction <= 1 ? 1 : fraction <= 2 ? 2 : fraction <= 5 ? 5 : 10)
    }
    readonly property real defaultDistance: sceneModel.frameRadius * 3.8
    readonly property real cameraDistance: camera.z
    readonly property vector3d cameraRotation: origin.eulerRotation
    readonly property vector3d cameraCenter: origin.scenePosition
    property string selectedTarget: "sun"
    readonly property vector3d targetPosition: {
        if (selectedTarget === "object") return sceneModel.objectPosition
        if (selectedTarget === "perihelion") return targetOrbit.perihelionPosition
        const index = sceneModel.planets.findIndex(planet => planet.name === selectedTarget)
        return index >= 0 && index < sceneModel.planetPositions.length ? sceneModel.planetPositions[index] : Qt.vector3d(0, 0, 0)
    }
    readonly property var legendEntries: {
        const entries = [{key: "sun", label: qsTr("Sun"), color: AppTheme.amber}]
        for (const planet of sceneModel.planets) entries.push({key: planet.name, label: planet.name, color: planetColor(planet.name)})
        entries.push({key: "object", label: solution.objectTitle, color: AppTheme.accent})
        if (showPerihelion) entries.push({key: "perihelion", label: qsTr("Perihelion"), color: AppTheme.red})
        return entries
    }
    readonly property var bodyLabels: {
        const entries = [
            {name: qsTr("Sun"), position: Qt.vector3d(0, 0, 0), color: AppTheme.amber, enabled: true},
            {name: solution.objectTitle, position: sceneModel.objectPosition, color: AppTheme.accent, enabled: true}
        ]
        for (let i = 0; i < sceneModel.planets.length; ++i) {
            const planet = sceneModel.planets[i]
            const available = i < sceneModel.planetPositions.length
            entries.push({name: planet.name, position: available ? sceneModel.planetPositions[i] : Qt.vector3d(0, 0, 0),
                          color: planetColor(planet.name), enabled: available})
        }
        return entries
    }
    color: AppTheme.nav
    clip: true

    function zoom(factor) { zoomAnimation.stop(); camera.z = Math.max(sceneModel.frameRadius * 0.15, Math.min(sceneModel.frameRadius * 80, camera.z * factor)) }
    function stopTransition() { centerAnimation.stop(); zoomAnimation.stop() }
    function resetView() { stopTransition(); selectedTarget = "sun"; origin.position = Qt.vector3d(0, 0, 0); origin.eulerRotation = Qt.vector3d(-52, 32, 0); camera.z = defaultDistance }
    function topView() { stopTransition(); origin.position = Qt.vector3d(0, 0, 0); origin.eulerRotation = Qt.vector3d(-90, 0, 0) }
    function follow(key) {
        stopTransition()
        const previousCenter = origin.scenePosition
        selectedTarget = key
        origin.position = previousCenter.minus(targetPosition)
        centerAnimation.restart()
        zoomAnimation.to = Math.max(sceneModel.frameRadius * 0.15, 1)
        zoomAnimation.restart()
    }
    function planetColor(name: string): color {
        if (name === "Earth") return AppTheme.blue
        if (name === "Mars") return AppTheme.red
        if (name === "Venus") return AppTheme.amber
        const planet = sceneModel.planets.find(planet => planet.name === name)
        return planet ? planet.color : AppTheme.muted
    }
    onShowPerihelionChanged: { if (!showPerihelion && selectedTarget === "perihelion") follow("object") }
    Connections { target: scene.sceneModel; function onChanged() { scene.resetView() } }
    Component.onCompleted: resetView()

    Vector3dAnimation { id: centerAnimation; target: origin; property: "position"; to: Qt.vector3d(0, 0, 0); duration: 500; easing.type: Easing.InOutQuad }
    NumberAnimation { id: zoomAnimation; target: camera; property: "z"; duration: 500; easing.type: Easing.InOutQuad }
    OrbitGrid {
        anchors.fill: parent
        visible: scene.showGrid
        sceneCamera: camera
        step: scene.gridStep
        center: scene.cameraCenter
    }
    View3D {
        id: view
        objectName: "orbit-viewport"
        anchors.fill: parent
        camera: camera
        environment: SceneEnvironment {
            backgroundMode: SceneEnvironment.Transparent
            antialiasingMode: SceneEnvironment.MSAA
            antialiasingQuality: SceneEnvironment.High
            tonemapMode: SceneEnvironment.TonemapModeLinear
        }
        Node {
            position: scene.targetPosition
            Node {
                id: origin
                objectName: "orbit-camera-origin"
                PerspectiveCamera {
                    id: camera
                    fieldOfView: 38
                    clipNear: Math.min(0.001, scene.sceneModel.frameRadius * 0.001)
                    clipFar: Math.max(100, scene.sceneModel.frameRadius * 150)
                    onZChanged: {
                        const bounded = Math.max(scene.sceneModel.frameRadius * 0.15, Math.min(scene.sceneModel.frameRadius * 80, z))
                        if (bounded !== z) z = bounded
                    }
                }
            }
        }
        Repeater3D {
            model: scene.sceneModel.planets
            Node {
                id: planet
                required property var modelData
                required property int index
                OrbitPath {
                    objectName: "orbit-path-" + planet.modelData.name
                    visible: scene.showPlanets
                    a: planet.modelData.a
                    e: planet.modelData.e
                    inclination: planet.modelData.i
                    node: planet.modelData.Om
                    peri: planet.modelData.w
                    color: UI.mix(scene.planetColor(planet.modelData.name), AppTheme.nav, 0.55)
                }
                OrbitBody {
                    objectName: "orbit-body-" + planet.modelData.name
                    visible: planet.index < scene.sceneModel.planetPositions.length
                    position: visible ? scene.sceneModel.planetPositions[planet.index] : Qt.vector3d(0, 0, 0)
                    radius: camera.z * 0.0024
                    color: scene.planetColor(planet.modelData.name)
                }
            }
        }
        OrbitPath {
            id: targetOrbit
            a: scene.solution.orbit.semiMajorAxis
            e: scene.solution.orbit.eccentricity
            q: scene.solution.orbit.perigee
            inclination: scene.solution.orbit.inclination
            node: scene.solution.orbit.ascendingNode
            peri: scene.solution.orbit.argPericentre
        }
        OrbitBody { color: AppTheme.amber; radius: camera.z * 0.0045 }
        OrbitBody { objectName: "orbit-object"; position: scene.sceneModel.objectPosition; radius: camera.z * 0.0028 }
        OrbitBody {
            visible: scene.showPerihelion
            position: targetOrbit.perihelionPosition
            radius: camera.z * 0.0014
            color: AppTheme.red
        }
    }
    OrbitCameraController {
        id: controller
        anchors.fill: parent
        origin: origin
        camera: camera
        automaticClipping: false
        mouseEnabled: !centerAnimation.running && !zoomAnimation.running
        acceptedButtons: Qt.LeftButton
        HoverHandler { cursorShape: controller.inputsNeedProcessing ? Qt.ClosedHandCursor : Qt.OpenHandCursor }
    }
    Repeater {
        model: scene.bodyLabels
        CustomText {
            id: label
            required property var modelData
            required property int index
            objectName: "orbit-label-" + modelData.name
            readonly property vector3d projected: {
                const cameraTransform = camera.sceneTransform
                const viewportWidth = view.width
                const viewportHeight = view.height
                return view.mapFrom3DScene(modelData.position)
            }
            readonly property bool collides: {
                const cameraTransform = camera.sceneTransform
                for (let i = 0; i < index; ++i) {
                    const other = scene.bodyLabels[i]
                    if (!other.enabled) continue
                    const point = view.mapFrom3DScene(other.position)
                    if (point.z > 0 && Math.abs(projected.x - point.x) < Math.max(width, other.name.length * 7) + 4
                        && Math.abs(projected.y - point.y) < 14) return true
                }
                return false
            }
            x: projected.x + 12
            y: projected.y - 18
            visible: scene.showLabels && modelData.enabled && !collides && projected.z > 0 && x > 8 && x + width < scene.width - 8 && y > 140 && y + height < legendFooter.y - 8
            width: Math.min(implicitWidth, 160)
            text: modelData.name
            color: modelData.color
            font.pixelSize: 11
        }
    }
    RowLayout {
        anchors.left: parent.left
        anchors.top: parent.top
        anchors.margins: 20
        spacing: 4
        CustomButton { objectName: "orbit-zoom-in"; iconName: "plus"; accessibleLabel: qsTr("Zoom in"); onClicked: scene.zoom(0.8) }
        CustomButton { objectName: "orbit-zoom-out"; iconName: "minus"; accessibleLabel: qsTr("Zoom out"); onClicked: scene.zoom(1.25) }
        CustomButton { objectName: "orbit-reset-view"; iconName: "rotate-ccw"; accessibleLabel: qsTr("Reset view"); onClicked: scene.resetView() }
        CustomButton { objectName: "orbit-top-view"; iconName: "view"; accessibleLabel: qsTr("Ecliptic view"); onClicked: scene.topView() }
    }
    Column {
        anchors.right: parent.right
        anchors.top: parent.top
        anchors.margins: 24
        width: 190
        spacing: 15
        CustomText { width: parent.width; text: scene.solution.objectTitle; horizontalAlignment: Text.AlignRight; font.weight: Font.DemiBold }
        KeyValueList {
            width: parent.width
            spacing: 8
            rows: [
                ["a", UI.number(scene.solution.orbit.semiMajorAxis, 6) + " AU"],
                ["e", UI.number(scene.solution.orbit.eccentricity, 6)],
                ["i", UI.number(scene.solution.orbit.inclination, 4) + "\u00b0"],
                [qsTr("Period"), UI.number(scene.solution.orbit.orbitalPeriod, 2) + " d"]
            ]
        }
    }
    ColumnLayout {
        id: legendFooter
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.bottom: parent.bottom
        anchors.margins: 20
        spacing: 8
        OrbitLegend {
            Layout.fillWidth: true
            entries: scene.legendEntries
            selectedTarget: scene.selectedTarget
            onTargetSelected: key => scene.follow(key)
        }
        RowLayout {
            Layout.fillWidth: true
            spacing: 16
            CustomText { text: scene.showGrid ? qsTr("Grid: %1 AU").arg(UI.number(scene.gridStep, scene.gridStep < 1 ? 1 : 0)) : ""; secondary: true; font.pixelSize: 11 }
            CustomText { Layout.fillWidth: true; text: scene.sceneModel.error; horizontalAlignment: Text.AlignRight; color: AppTheme.amber; font.pixelSize: 11; wrapMode: Text.Wrap }
        }
    }
}
