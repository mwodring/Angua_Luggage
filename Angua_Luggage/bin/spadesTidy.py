from ..LuggageInterface import spadesTidy
import sys, logging, logging.config, os

logging.basicConfig(stream = sys.stdout, level=logging.INFO)
LOG = logging.getLogger(__name__)

def main():
	in_dir = os.path.abspath(sys.argv[1])
	out_dir = os.path.abspath(sys.argv[2])
	tidier = spadesTidy("in", in_dir)
	tidier.spadesToDir(out_dir)
	
if __name__ == "__main__":
    sys.exit(main())