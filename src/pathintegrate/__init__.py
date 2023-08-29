from pkg_resources import get_distribution
__version__ = get_distribution('pathintegrate').version

from .pathintegrate import PathIntegrate
from .app import launch_network_app