import numpy as np
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
                              QLabel, QRadioButton, QButtonGroup, QComboBox,
                              QFrame)
from PyQt5.QtCore import pyqtSignal

_STYLE_ACTIVE   = "QPushButton { background-color: #2a5db0; color: white; border-radius: 4px; padding: 4px 8px; }"
_STYLE_INACTIVE = ""


class MeasurementWidget(QWidget):
    measure_select_toggled = pyqtSignal(bool)
    frame_type_changed     = pyqtSignal()
    set_endpoint_active    = pyqtSignal(int)
    clear_requested        = pyqtSignal()
    axis_frame_changed     = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.frame_types = ['center', 'center']
        self.axis_frame  = 'global'
        self._current_ep = -1  # -1=regular, 0=select A, 1=select B

        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(4, 4, 4, 4)
        main_layout.setSpacing(4)

        # Three-column layout: Point A | Point B | Results
        columns_layout = QHBoxLayout()
        columns_layout.setSpacing(0)

        self.point_label: list[QLabel]       = []
        self.coord_label: list[QLabel]       = []
        self.frame_group: list[QButtonGroup] = []
        self.set_btn:     list[QPushButton]  = []

        endpoint_names = ["Point A", "Point B"]
        for ep in range(2):
            if ep > 0:
                columns_layout.addWidget(self._vline())

            col_widget = QWidget()
            col_layout = QVBoxLayout(col_widget)
            col_layout.setContentsMargins(6, 4, 6, 4)
            col_layout.setSpacing(3)

            btn = QPushButton(f"Select {endpoint_names[ep]}")
            btn.clicked.connect(lambda _, e=ep: self._on_select_btn_clicked(e))
            btn.setStyleSheet(_STYLE_INACTIVE)
            self.set_btn.append(btn)
            col_layout.addWidget(btn)

            point_lbl = QLabel("Click in viewport…")
            point_lbl.setWordWrap(True)
            self.point_label.append(point_lbl)
            col_layout.addWidget(point_lbl)

            coord_lbl = QLabel("")
            self.coord_label.append(coord_lbl)
            col_layout.addWidget(coord_lbl)

            radio_row = QHBoxLayout()
            grp = QButtonGroup(self)
            self.frame_group.append(grp)
            for i, name in enumerate(['Back', 'Center', 'Front']):
                rb = QRadioButton(name)
                if name == 'Center':
                    rb.setChecked(True)
                grp.addButton(rb, i)
                radio_row.addWidget(rb)
            grp.idClicked.connect(
                lambda btn_id, e=ep: self._on_frame_radio_clicked(e, btn_id)
            )
            col_layout.addLayout(radio_row)

            col_layout.addStretch()
            columns_layout.addWidget(col_widget)

        # Results column
        columns_layout.addWidget(self._vline())

        results_widget = QWidget()
        results_layout = QVBoxLayout(results_widget)
        results_layout.setContentsMargins(6, 4, 6, 4)
        results_layout.setSpacing(4)

        self.clear_btn = QPushButton("Clear Measurement")
        self.clear_btn.clicked.connect(self.clear_requested.emit)
        results_layout.addWidget(self.clear_btn)

        self.distance_label = QLabel("Distance: —")
        results_layout.addWidget(self.distance_label)

        axis_row = QHBoxLayout()
        axis_row.addWidget(QLabel("Axis frame:"))
        self.axis_combo = QComboBox()
        self.axis_combo.addItems(["Global", "Local A", "Local B"])
        self.axis_combo.currentIndexChanged.connect(self._on_axis_combo_changed)
        axis_row.addWidget(self.axis_combo)
        results_layout.addLayout(axis_row)

        self.dx_label = QLabel("dX: —")
        self.dy_label = QLabel("dY: —")
        self.dz_label = QLabel("dZ: —")
        results_layout.addWidget(self.dx_label)
        results_layout.addWidget(self.dy_label)
        results_layout.addWidget(self.dz_label)
        results_layout.addStretch()

        columns_layout.addWidget(results_widget)
        main_layout.addLayout(columns_layout)

    def _vline(self):
        line = QFrame()
        line.setFrameShape(QFrame.VLine)
        line.setStyleSheet("color: #888;")
        return line

    def _set_mode(self, ep: int):
        self._current_ep = ep
        for i in range(2):
            self.set_btn[i].setStyleSheet(_STYLE_ACTIVE if i == ep else _STYLE_INACTIVE)

    def _on_select_btn_clicked(self, ep: int):
        if self._current_ep == ep:
            self._set_mode(-1)
            self.measure_select_toggled.emit(False)
        else:
            was_inactive = (self._current_ep == -1)
            self._set_mode(ep)
            if was_inactive:
                self.measure_select_toggled.emit(True)
            self.set_endpoint_active.emit(ep)

    def deactivate(self):
        self._set_mode(-1)

    def _on_frame_radio_clicked(self, ep, btn_id):
        types = ['proximal', 'center', 'distal']
        self.frame_types[ep] = types[btn_id]
        self.frame_type_changed.emit()

    def _on_axis_combo_changed(self, index):
        options = ['global', 'local_a', 'local_b']
        self.axis_frame = options[index]
        self.axis_frame_changed.emit()

    def set_active_endpoint(self, ep: int):
        self._set_mode(ep)

    def set_point_info(self, ep: int, joint_index: int, type_name: str, pos):
        self.point_label[ep].setText(f"Joint {joint_index} ({type_name})")
        if pos is not None:
            self.coord_label[ep].setText(f"({pos[0]:.1f}, {pos[1]:.1f}, {pos[2]:.1f})")
        else:
            self.coord_label[ep].setText("")

    def set_result(self, dist, components):
        if dist is None:
            self.distance_label.setText("Distance: —")
            self.dx_label.setText("dX: —")
            self.dy_label.setText("dY: —")
            self.dz_label.setText("dZ: —")
        else:
            self.distance_label.setText(f"Distance: {dist:.2f} mm")
            self.dx_label.setText(f"dX: {components[0]:.2f}")
            self.dy_label.setText(f"dY: {components[1]:.2f}")
            self.dz_label.setText(f"dZ: {components[2]:.2f}")

    def reset(self):
        for ep in range(2):
            self.point_label[ep].setText("Click in viewport…")
            self.coord_label[ep].setText("")
            center_btn = self.frame_group[ep].button(0)
            if center_btn is not None:
                center_btn.setChecked(True)
            self.frame_types[ep] = 'center'
        self.set_result(None, None)
        self._set_mode(-1)
