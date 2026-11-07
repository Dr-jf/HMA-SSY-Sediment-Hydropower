"""Run historical annual-capacity updates and signed paired generation loss."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from hma_models.hydropower import sediment_cli


if __name__ == "__main__":
    sediment_cli()
