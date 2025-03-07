# -*- mode: python ; coding: utf-8 -*-
block_cipher = None
a = Analysis(
    ['GUI.py'],
    pathex=['.'],
    binaries=[],
    datas=[
        ('ui/rotate_unchecked.png', 'ui'),
        ('ui/rotate_checked.png', 'ui'),
        ('ui/move_unchecked.png', 'ui'),
        ('ui/move_checked.png', 'ui'),
        ('save/test.tree', 'save'),
        ('save/test.chain', 'save'),
        ('save/r.dxf', 'save'),
        ('save/r.chain', 'save'),
        ('save/autosave/autosave_0.chain', 'save/autosave'),
        ('referenceMeshes/stanfordBunnyLowPoly.stl', 'referenceMeshes'),
        ('referenceMeshes/meshSources.txt', 'referenceMeshes'),
        ('referenceMeshes/legBones.stl', 'referenceMeshes'),
        ('referenceMeshes/humanArmBones.stl', 'referenceMeshes'),
        ('referenceMeshes/arm.stl', 'referenceMeshes'),
    ],
    hiddenimports=[],
    hookspath=[],
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)
exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='GUI',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False
)
coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    name='GUI'
)