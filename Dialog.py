"""
@author: Raymond Feng, Andy Wang, Daniel Feshbach
"""

from PyQt5.QtWidgets import QVBoxLayout, QPushButton, QLabel, QDialog, QLineEdit


class EditJointStateDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Edit Joint State')
        self.setGeometry(100, 100, 200, 100)

        layout = QVBoxLayout()

        self.state_input = QLineEdit(self)
        self.state_input.setPlaceholderText('Enter new joint state')
        layout.addWidget(QLabel('State:'))
        layout.addWidget(self.state_input)

        self.apply_button = QPushButton('Apply', self)
        self.apply_button.clicked.connect(self.onApplyClicked)
        layout.addWidget(self.apply_button)

        self.cancel_button = QPushButton('Cancel', self)
        self.cancel_button.clicked.connect(self.reject)
        layout.addWidget(self.cancel_button)

        self.setLayout(layout)

    def onApplyClicked(self):
        self.accept()

    def get_state(self):
        try:
            state = float(self.state_input.text())
            return state
        except ValueError:
            self.show_error("Please enter a valid join state.")
            # QMessageBox.warning(self, "Invalid Input", "Please enter a valid join state.")
            self.exec_() 
            return None