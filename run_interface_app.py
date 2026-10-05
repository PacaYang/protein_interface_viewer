#!/usr/bin/env python3
"""Start the local protein-interface analysis app."""

import argparse
import os

import uvicorn


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--data-dir", default=".interface_app_data")
    args = parser.parse_args()
    os.environ["INTERFACE_APP_DATA"] = args.data_dir
    uvicorn.run("interface_app.app:app", host=args.host, port=args.port, reload=False)


if __name__ == "__main__":
    main()
