"""Run historical baseline hydropower from prepared private CSV inputs."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from hma_models.hydropower import baseline_cli


if __name__ == "__main__":
    baseline_cli()
