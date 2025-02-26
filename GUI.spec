<<<<<<< Updated upstream
# my_project.spec
# PyInstaller spec file example

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
=======
# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['GUI.py'],
    pathex=[],
    binaries=[],
    datas=[],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)
>>>>>>> Stashed changes

exe = EXE(
    pyz,
    a.scripts,
<<<<<<< Updated upstream
    [],
    exclude_binaries=True,
    name='GUI',  
=======
    a.binaries,
    a.datas,
    [],
    name='GUI',
>>>>>>> Stashed changes
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
<<<<<<< Updated upstream
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
=======
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
>>>>>>> Stashed changes
)
