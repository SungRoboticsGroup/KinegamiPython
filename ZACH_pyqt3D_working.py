import sys
from PyQt5.QtCore import Qt, QSize
from PyQt5.QtGui import QColor, QVector3D, QQuaternion
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QHBoxLayout
from PyQt5.Qt3DCore import QEntity, QTransform
from PyQt5.Qt3DExtras import Qt3DWindow, QCuboidMesh, QPhongMaterial, QOrbitCameraController, QConeMesh, QCylinderMesh
from PyQt5.Qt3DRender import QPointLight, QObjectPicker, QRayCaster

class TranslationGizmo(QEntity):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.arrows = {}
        self.pickers = {}
        self.transforms = {}
        self.create_arrows()

    def create_arrows(self):
        axes = {
            'x': {'color': QColor(255, 0, 0), 'rotation': QQuaternion.fromAxisAndAngle(QVector3D(0, 0, 1), -90)},
            'y': {'color': QColor(0, 255, 0), 'rotation': QQuaternion()},
            'z': {'color': QColor(0, 0, 255), 'rotation': QQuaternion.fromAxisAndAngle(QVector3D(1, 0, 0), 90)},
        }

        for axis, props in axes.items():
            arrow_entity = QEntity(self)

            cylinder_mesh = QCylinderMesh()
            cylinder_mesh.setRadius(0.1)
            cylinder_mesh.setLength(0.75)
            cylinder_mesh.setRings(10)
            cylinder_mesh.setSlices(20)

            transform = QTransform()
            transform.setRotation(props['rotation'])
            transform.setTranslation(QVector3D(*{'x': (1, 0, 0), 'y': (0, 1, 0), 'z': (0, 0, 1)}[axis]))

            material = QPhongMaterial()
            material.setDiffuse(props['color'])

            picker = QObjectPicker(arrow_entity)
            picker.setHoverEnabled(True)
            picker.setDragEnabled(True)

            arrow_entity.addComponent(cylinder_mesh)
            arrow_entity.addComponent(transform)
            arrow_entity.addComponent(material)
            arrow_entity.addComponent(picker)

            self.arrows[axis] = arrow_entity
            self.pickers[axis] = picker
            self.transforms[axis] = transform

    def setPosition(self, position):
        for axis, transform in self.transforms.items():
            offset = QVector3D(*{'x': (1, 0, 0), 'y': (0, 1, 0), 'z': (0, 0, 1)}[axis])
            transform.setTranslation(position + offset)

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Translation Gizmo Demo")
        self.resize(800, 600)

        self.view = Qt3DWindow()
        self.container = self.createWindowContainer(self.view)
        self.container.setMinimumSize(QSize(800, 600))
        self.container.setFocusPolicy(Qt.TabFocus)

        self.widget = QWidget()
        self.layout = QHBoxLayout(self.widget)
        self.layout.addWidget(self.container)
        self.setCentralWidget(self.widget)

        self.root_entity = QEntity()
        self.view.setRootEntity(self.root_entity)

        self.camera = self.view.camera()
        self.camera.lens().setPerspectiveProjection(45.0, self.width() / self.height(), 0.1, 1000)
        self.camera.setPosition(QVector3D(4, 4, 4))
        self.camera.setViewCenter(QVector3D(0, 0, 0))

        # self.cam_controller = QOrbitCameraController(self.root_entity)
        # self.cam_controller.setLinearSpeed(50)
        # self.cam_controller.setLookSpeed(180)
        # self.cam_controller.setCamera(self.camera)

        self.light_entity = QEntity(self.root_entity)
        self.light = QPointLight(self.light_entity)
        self.light.setColor(QColor(255, 255, 255))
        self.light.setIntensity(1)
        self.light_transform = QTransform()
        self.light_transform.setTranslation(QVector3D(10, 10, 10))
        self.light_entity.addComponent(self.light)
        self.light_entity.addComponent(self.light_transform)

        self.cube_entity = QEntity(self.root_entity)
        self.cube_mesh = QCuboidMesh()
        self.cube_transform = QTransform()
        self.cube_material = QPhongMaterial()
        self.cube_material.setDiffuse(QColor(200, 200, 200))
        self.cube_entity.addComponent(self.cube_mesh)
        self.cube_entity.addComponent(self.cube_transform)
        self.cube_entity.addComponent(self.cube_material)

        self.gizmo = TranslationGizmo(self.root_entity)
        self.selected_axis = None
        self.start_pos_3D = None

        for axis, picker in self.gizmo.pickers.items():
            picker.pressed.connect(lambda event, a=axis: self.on_picker_pressed(event, a))
            picker.moved.connect(self.on_picker_moved)
            picker.released.connect(self.on_picker_released)

    def get_world_coordinates(self, event):
        pos = event.position()
        ndc_x = (2.0 * pos.x()) / self.view.width() - 1.0
        ndc_y = 1.0 - (2.0 * pos.y()) / self.view.height()
        ndc = QVector3D(ndc_x, ndc_y, -1.0)

        inverted_matrix = (self.camera.projectionMatrix() * self.camera.viewMatrix()).inverted()[0]
        near_point = inverted_matrix.map(ndc)
        ndc.setZ(1.0)
        far_point = inverted_matrix.map(ndc)

        direction = far_point - near_point
        direction.normalize()

        return near_point, direction

    def on_picker_pressed(self, event, axis):
        self.selected_axis = axis

        origin, direction = self.get_world_coordinates(event)
        self.start_pos_3D = self.compute_intersection(origin, direction, axis)

        event.setAccepted(True)

    def on_picker_moved(self, event):
        if self.selected_axis:
            origin, direction = self.get_world_coordinates(event)

            new_pos_3D = self.compute_intersection(origin, direction, self.selected_axis)

            if new_pos_3D and self.start_pos_3D:
                delta_vector = new_pos_3D - self.start_pos_3D
                axis_vector = self.get_axis_vector(self.selected_axis)
                movement = QVector3D.dotProduct(delta_vector, axis_vector.normalized())

                translation = self.cube_transform.translation()
                translation += axis_vector.normalized() * movement
                self.cube_transform.setTranslation(translation)
                self.gizmo.setPosition(translation)

                self.start_pos_3D = new_pos_3D

            event.setAccepted(True)

    def on_picker_released(self, event):
        self.selected_axis = None
        self.start_pos_3D = None
        event.setAccepted(True)

    def compute_intersection(self, origin, direction, axis):
        cube_pos = self.cube_transform.translation()
        axis_vector = self.get_axis_vector(axis)
        plane_normal = QVector3D.crossProduct(
            QVector3D.crossProduct(direction, axis_vector), axis_vector
        ).normalized()

        denom = QVector3D.dotProduct(direction, plane_normal)
        if abs(denom) > 1e-6:
            t = QVector3D.dotProduct(cube_pos - origin, plane_normal) / denom
            intersection_point = origin + direction * t
            return intersection_point
        else:
            return None

    def get_axis_vector(self, axis):
        if axis == 'x':
            return QVector3D(1, 0, 0)
        elif axis == 'y':
            return QVector3D(0, 1, 0)
        elif axis == 'z':
            return QVector3D(0, 0, 1)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())