import pyqtgraph as pg
import pyqtgraph.opengl as gl
import numpy as np

class MyGLV(gl.GLViewWidget):
    def mouseReleaseEvent(self, ev):
        lpos = ev.position() if hasattr(ev, 'position') else ev.localPos()
        region = [lpos.x()-5, lpos.y()-5, 10, 10]
        # itemsAt seems to take in device pixels
        dpr = self.devicePixelRatioF()
        region = tuple([x * dpr for x in region])
        for item in self.itemsAt(region):
            print(item.objectName())

pg.mkQApp()

glv = MyGLV()
# glv.setCameraParams(elevation=90, azimuth=-90, distance=50)
# X points right, Y points up
glv.show()

side = 4
names = ['red', 'green', 'blue']
for idx, name in enumerate(names):
    box = np.zeros((side, side, 4), dtype=np.uint8)
    box[..., [idx, 3]] = 255
    img = gl.GLImageItem(box)
    img.setObjectName(name)
    img.translate(-side/2, -side/2, 0)  # center the box
    img.translate((idx-1)*side*2, (idx-1)*side*1, 0)
    glv.addItem(img)

pg.exec()