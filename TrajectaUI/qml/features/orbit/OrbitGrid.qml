import QtQuick
import QtQuick3D
import TrajectaUI
import "../../shared/UI.js" as UI

View3D {
    id: grid
    objectName: "orbit-grid-viewport"
    required property PerspectiveCamera sceneCamera
    required property real step
    required property vector3d center
    readonly property real visibleHalfWidth: sceneCamera.z * Math.tan(sceneCamera.fieldOfView * Math.PI / 360)
        * Math.min(8, Math.max(1, width / Math.max(1, height)))
    readonly property int halfDivisions: Math.max(8, Math.ceil(visibleHalfWidth / step) + 2)
    camera: gridCamera
    environment: SceneEnvironment {
        backgroundMode: SceneEnvironment.Color
        clearColor: AppTheme.nav
        antialiasingMode: SceneEnvironment.MSAA
        antialiasingQuality: SceneEnvironment.High
        tonemapMode: SceneEnvironment.TonemapModeLinear
    }
    PerspectiveCamera {
        id: gridCamera
        position: grid.sceneCamera.scenePosition
        rotation: grid.sceneCamera.sceneRotation
        fieldOfView: grid.sceneCamera.fieldOfView
        clipNear: grid.sceneCamera.clipNear
        clipFar: grid.sceneCamera.clipFar
    }
    Model {
        objectName: "orbit-grid"
        position: Qt.vector3d(Math.round(grid.center.x / grid.step) * grid.step, 0,
                              Math.round(grid.center.z / grid.step) * grid.step)
        geometry: GridGeometry { size: grid.step * grid.halfDivisions * 2; divisions: grid.halfDivisions * 2 }
        materials: PrincipledMaterial {
            lighting: PrincipledMaterial.NoLighting
            baseColor: UI.mix(AppTheme.line, AppTheme.nav, 0.45)
        }
        castsShadows: false
        receivesShadows: false
    }
}
