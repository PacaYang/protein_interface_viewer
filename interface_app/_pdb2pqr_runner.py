"""Run managed PDB2PQR wheels using the app's Python and NumPy installation."""

import json
import sys


if __name__ == "__main__":
    sys.path[:0] = json.loads(sys.argv.pop(1))
    from pdb2pqr.main import main

    main()
