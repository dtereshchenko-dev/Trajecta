pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import QtQuick3D
import QtQuick3D.Helpers
import TrajectaUI
import "../../shared/UI.js" as UI

Rectangle {
    id: scene
    objectName: "relative-orbit-scene"
    required property SolutionViewModel solution
    readonly property OrbitSceneModel model: solution.scene
    readonly property RelativeMotionModel relative: model.relative
    readonly property real defaultDistance: Math.max(model.frameRadius * 3.8, 0.0001)
    readonly property color distantEnergyColor: AppTheme.name === "light" || AppTheme.name === "dim" ? "#715887" : "#b098c6"
    property bool showHill: true
    property bool showMoon: true
    property bool orthographicView: false
    readonly property real orthographicMagnification: UI.orthographicMagnification(view.height, camera.z, camera.fieldOfView)
    color: AppTheme.nav
    clip: true

    function resetView() {
        centerAnimation.stop()
        focusZoomAnimation.stop()
        orthographicView = false
        origin.position = Qt.vector3d(0, 0, 0)
        origin.eulerRotation = Qt.vector3d(-52, 32, 0)
        camera.z = defaultDistance
    }
    function topView() { centerAnimation.stop(); focusZoomAnimation.stop(); origin.position = Qt.vector3d(0, 0, 0); origin.eulerRotation = Qt.vector3d(-90, 0, 0); orthographicView = true }
    function focusOnObject() {
        centerAnimation.stop()
        focusZoomAnimation.stop()
        centerAnimation.to = scene.model.objectPosition
        focusZoomAnimation.to = Math.max(scene.relative.relativePlanetRadius * 20, scene.relative.relativeHillRadius * 2.5)
        centerAnimation.restart()
        focusZoomAnimation.restart()
    }
    function zoom(factor) { focusZoomAnimation.stop(); camera.z = Math.max(relative.relativePlanetRadius * 2.5, Math.min(model.frameRadius * 80, camera.z * factor)) }
    Component.onCompleted: resetView()
    Connections { target: scene.model; function onChanged() { scene.resetView() } }
    Vector3dAnimation { id: centerAnimation; target: origin; property: "position"; duration: 500; easing.type: Easing.InOutQuad }
    NumberAnimation { id: focusZoomAnimation; target: camera; property: "z"; duration: 500; easing.type: Easing.InOutQuad }

    OrbitGrid {
        anchors.fill: parent
        visible: scene.relative.relativeAvailable
        sceneCamera: camera
        step: Math.pow(10, Math.floor(Math.log10(Math.max(1e-8, camera.z / 6))))
        center: origin.scenePosition
        orthographic: scene.orthographicView
        orthographicMagnification: scene.orthographicMagnification
    }
    View3D {
        id: view
        objectName: "relative-viewport"
        anchors.fill: parent
        camera: scene.orthographicView ? orthographicCamera : camera
        environment: SceneEnvironment {
            backgroundMode: SceneEnvironment.Transparent
            antialiasingMode: SceneEnvironment.MSAA
            antialiasingQuality: SceneEnvironment.High
            tonemapMode: SceneEnvironment.TonemapModeLinear
        }
        Node {
            id: origin
            objectName: "relative-camera-origin"
            onEulerRotationChanged: {
                if (scene.orthographicView && (Math.abs(eulerRotation.x + 90) > 0.01 || Math.abs(eulerRotation.y) > 0.01 || Math.abs(eulerRotation.z) > 0.01))
                    scene.orthographicView = false
            }
            PerspectiveCamera {
                id: camera
                objectName: "relative-camera"
                fieldOfView: 38
                clipNear: Math.max(1e-12, camera.z * 0.00001)
                clipFar: Math.max(1, scene.model.frameRadius * 150)
            }
            OrthographicCamera {
                id: orthographicCamera
                objectName: "relative-orthographic-camera"
                z: camera.z
                clipNear: camera.clipNear
                clipFar: camera.clipFar
                horizontalMagnification: scene.orthographicMagnification
                verticalMagnification: scene.orthographicMagnification
            }
        }
        Model {
            visible: scene.relative.relativeAvailable
            geometry: PolylineGeometry { independentSegments: true; points: scene.relative.relativeUnboundSegments }
            materials: PrincipledMaterial { lighting: PrincipledMaterial.NoLighting; baseColor: AppTheme.accent }
        }
        Model {
            visible: scene.relative.relativeAvailable
            geometry: PolylineGeometry { independentSegments: true; points: scene.relative.relativeBoundSegments }
            materials: PrincipledMaterial { lighting: PrincipledMaterial.NoLighting; baseColor: AppTheme.amber }
        }
        Model {
            visible: scene.relative.relativeAvailable
            geometry: PolylineGeometry { independentSegments: true; points: scene.relative.relativeDistantNegativeSegments }
            materials: PrincipledMaterial { lighting: PrincipledMaterial.NoLighting; baseColor: scene.distantEnergyColor }
        }
        Model {
            visible: scene.showMoon && scene.relative.relativeMoonPath.length > 0
            geometry: PolylineGeometry { points: scene.relative.relativeMoonPath }
            materials: PrincipledMaterial { lighting: PrincipledMaterial.NoLighting; baseColor: AppTheme.muted }
        }
        OrbitBody { objectName: "relative-planet"; visible: scene.relative.relativeAvailable; radius: Math.max(scene.relative.relativePlanetRadius, camera.z * 0.004); color: AppTheme.blue }
        OrbitBody { visible: scene.relative.relativeAvailable; position: scene.model.objectPosition; radius: camera.z * 0.003; color: AppTheme.accent }
        OrbitBody { visible: scene.showMoon && scene.relative.relativeMoonPath.length > 0; position: scene.relative.relativeMoon; radius: camera.z * 0.002; color: AppTheme.muted }
        RelativeHillGuides { model: scene.relative; camera: camera; showHill: scene.showHill }
        Repeater3D {
            model: scene.relative.relativeCaptureMarkers
            OrbitBody {
                required property var modelData
                position: modelData.position
                radius: camera.z * 0.0025
                color: modelData.start ? AppTheme.blue : AppTheme.red
            }
        }
    }
    OrbitCameraController {
        anchors.fill: parent
        origin: origin
        camera: camera
        automaticClipping: false
        mouseEnabled: !centerAnimation.running && !focusZoomAnimation.running
        acceptedButtons: Qt.LeftButton
    }
    OrbitViewControls {
        anchors.left: parent.left
        anchors.top: parent.top
        anchors.margins: 20
        objectPrefix: "relative"
        onZoomIn: scene.zoom(0.8)
        onZoomOut: scene.zoom(1.25)
        onResetView: scene.resetView()
        onTopView: scene.topView()
    }
    RelativeSunIndicator {
        anchors.left: parent.left
        anchors.top: parent.top
        anchors.leftMargin: 20
        anchors.topMargin: 62
        model: scene.relative
        camera: camera
    }
    Repeater {
        model: ["L1", "L2"]
        CustomText {
            required property string modelData
            readonly property vector3d projected: {
                if (!scene.relative.relativeAvailable || !view.camera || view.width <= 0 || view.height <= 0)
                    return Qt.vector3d(0, 0, -1)
                const cameraTransform = view.camera.sceneTransform
                const magnification = scene.orthographicView ? orthographicCamera.verticalMagnification : camera.fieldOfView
                const viewportWidth = view.width
                const viewportHeight = view.height
                const sign = modelData === "L1" ? 1 : -1
                return view.mapFrom3DScene(scene.relative.relativeSunDirection.times(sign * scene.relative.relativeHillRadius))
            }
            x: projected.x + 8
            y: projected.y - 15
            text: modelData + " \u2248"
            color: AppTheme.amber
            font.pixelSize: 11
            visible: scene.showHill && scene.relative.relativeAvailable && projected.z > 0 &&
                     x >= 0 && x + width < scene.width && y >= 0 && y + height < scene.height
        }
    }
    Repeater {
        model: scene.relative.relativeCaptureMarkers
        CustomText {
            required property var modelData
            readonly property vector3d projected: {
                if (!view.camera || view.width <= 0 || view.height <= 0) return Qt.vector3d(0, 0, -1)
                const cameraTransform = view.camera.sceneTransform
                const magnification = scene.orthographicView ? orthographicCamera.verticalMagnification : camera.fieldOfView
                const viewportWidth = view.width
                const viewportHeight = view.height
                return view.mapFrom3DScene(modelData.position)
            }
            x: projected.x + 9
            y: projected.y - 17
            text: modelData.start ? qsTr("Capture start %1").arg(modelData.event) : qsTr("Capture end %1").arg(modelData.event)
            color: modelData.start ? AppTheme.blue : AppTheme.red
            font.pixelSize: 11
            visible: projected.z > 0 && x >= 0 && x + width < scene.width && y >= 0 && y + height < scene.height
        }
    }
    CustomText {
        anchors.centerIn: parent
        visible: !scene.relative.relativeAvailable && !scene.relative.relativeBusy
        text: qsTr("Choose a planet and interval, then calculate the trajectory.")
        secondary: true
    }
    CustomText {
        anchors.centerIn: parent
        visible: scene.relative.relativeBusy
        text: qsTr("Calculating relative motion...")
        secondary: true
    }
    Flow {
        objectName: "relative-trajectory-legend"
        anchors.left: parent.left
        anchors.bottom: parent.bottom
        anchors.margins: 20
        width: Math.max(0, parent.width - 40)
        spacing: 14
        visible: scene.relative.relativeAvailable
        CustomText { text: scene.relative.relativeCenter; color: AppTheme.blue; font.pixelSize: 11 }
        CustomText { text: qsTr("Energy \u2265 0"); color: AppTheme.accent; font.pixelSize: 11 }
        CustomText { text: qsTr("Energy < 0, r < 3 Hill radii"); color: AppTheme.amber; font.pixelSize: 11 }
        CustomText { text: qsTr("Energy < 0, r \u2265 3 Hill radii"); color: scene.distantEnergyColor; font.pixelSize: 11 }
        CustomText { visible: scene.showMoon && scene.relative.relativeMoonPath.length > 0; text: qsTr("Moon"); secondary: true; font.pixelSize: 11 }
        CustomText { visible: scene.showHill; text: qsTr("Hill radius (approx.)"); secondary: true; font.pixelSize: 11 }
    }
}
