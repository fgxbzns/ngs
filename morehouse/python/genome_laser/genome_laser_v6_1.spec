# -*- mode: python -*-
a = Analysis(['genome_laser_v6_1.py'],
             pathex=['/home/guoxing/disk2/ngs/morehouse/python/genome_laser'],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)
pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='genome_laser_v6_1',
          debug=False,
          strip=None,
          upx=True,
          console=True )
