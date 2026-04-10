"""
Use pyAvaCore as a Python library to download EAWS bulletins.
Test with a single recent date to understand the data structure.
"""
import sys
sys.path.insert(0, '.')

# Import pyAvaCore
try:
    from avacore import pyAvaCore
    print("Imported pyAvaCore successfully")
except ImportError as e:
    print("Import error: {}".format(e))
    # Try alternate import
    try:
        import avacore
        print("avacore dir: {}".format(dir(avacore)))
    except:
        pass

# Check available modules/classes
try:
    from avacore import pyAvaCore
    print("pyAvaCore dir: {}".format([x for x in dir(pyAvaCore) if not x.startswith('_')]))
except:
    pass

# Try to download bulletins for a single day
try:
    from avacore.pyAvaCore import get_reports
    print("\nget_reports function found")
    
    # Get bulletins for 2024-01-15
    reports = get_reports(date='2024-01-15', regions=['AT-07'])
    print("Reports type: {}".format(type(reports)))
    if reports:
        print("Got {} reports".format(len(reports)))
        for r in reports[:3]:
            print("  Type: {}, attrs: {}".format(type(r).__name__, dir(r)[:10]))
except Exception as e:
    print("get_reports error: {}".format(e))

# Alternative: try the processor directly
try:
    from avacore.processor import process
    print("\nprocess function found")
except ImportError:
    pass

# Try the main entry point
try:
    from avacore.__main__ import main
    print("\nmain function found")
except ImportError:
    pass

# List all avacore submodules
import pkgutil
import avacore
print("\navacore submodules:")
for importer, modname, ispkg in pkgutil.walk_packages(avacore.__path__, avacore.__name__+'.'):
    print("  {} (pkg={})".format(modname, ispkg))
