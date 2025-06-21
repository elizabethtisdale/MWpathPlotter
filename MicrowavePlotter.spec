<<<<<<< HEAD

# -*- mode: python ; coding: utf-8 -*-

block_cipher = None

a = Analysis(['MWPathPainter.py'],
             pathex=[],
             
             binaries=[
                ('C:\\Users\\OEE2024_05\\anaconda3\\envs\\plotter_env\\Library\\bin\\gdal.dll', '.'),
                ('C:\\Users\\OEE2024_05\\anaconda3\\envs\\plotter_env\\Library\\bin\\proj_9.dll', '.'),
                ('C:\\Users\\OEE2024_05\\anaconda3\\envs\\plotter_env\\Library\\bin\\geos_c.dll', '.')
             ],
            
             datas=[
                ('C:\\Users\\OEE2024_05\\anaconda3\\envs\\plotter_env\\Library\\share\\proj', 'proj'),
                ('C:\\Users\\OEE2024_05\\anaconda3\\envs\\plotter_env\\Library\\share\\gdal', 'gdal')
             ],
             
             
             hiddenimports=[
                'openpyxl.cell._writer',
                'pyproj',
                'geopandas',
                'fiona',
                'pyogrio',
                'shapely',
                'pandas',
                'matplotlib',
                'tkinter',
				'fiona',
                'fiona.schema',
                'fiona.hooks',
                'fiona.drvsupport',
                'fiona.collection',
                'pyogrio',
                'pyogrio._geometry',
                'pyogrio._io',
                'pyogrio.errors'
             ],
             hookspath=[],
             hooksconfig={},
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)

pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)


exe = EXE(pyz,
          a.scripts,
          a.binaries, 
          a.zipfiles, 
          a.datas,    
          [],
          name='MWpathPainter',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          runtime_tmpdir=None,
          
          console=False,
          # icon='for a future icon'
          )
=======

# -*- mode: python ; coding: utf-8 -*-

block_cipher = None

a = Analysis(['MWPathPainter.py'],
             pathex=[],
             
             binaries=[
                ('C:\\Users\\OEE2024_05\\anaconda3\\envs\\plotter_env\\Library\\bin\\gdal.dll', '.'),
                ('C:\\Users\\OEE2024_05\\anaconda3\\envs\\plotter_env\\Library\\bin\\proj_9.dll', '.'),
                ('C:\\Users\\OEE2024_05\\anaconda3\\envs\\plotter_env\\Library\\bin\\geos_c.dll', '.')
             ],
            
             datas=[
                ('C:\\Users\\OEE2024_05\\anaconda3\\envs\\plotter_env\\Library\\share\\proj', 'proj'),
                ('C:\\Users\\OEE2024_05\\anaconda3\\envs\\plotter_env\\Library\\share\\gdal', 'gdal')
             ],
             
             
             hiddenimports=[
                'openpyxl.cell._writer',
                'pyproj',
                'geopandas',
                'fiona',
                'pyogrio',
                'shapely',
                'pandas',
                'matplotlib',
                'tkinter',
				'fiona',
                'fiona.schema',
                'fiona.hooks',
                'fiona.drvsupport',
                'fiona.collection',
                'pyogrio',
                'pyogrio._geometry',
                'pyogrio._io',
                'pyogrio.errors'
             ],
             hookspath=[],
             hooksconfig={},
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)

pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)


exe = EXE(pyz,
          a.scripts,
          a.binaries, 
          a.zipfiles, 
          a.datas,    
          [],
          name='MWpathPainter',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          runtime_tmpdir=None,
          
          console=False,
          # icon='for a future icon'
          )
>>>>>>> 3683fbb36275adddd3eb8cb0935a342a7d7dd4d2
