try:
    from .fortran_backend import build_extension
except ImportError:
    from fortran_backend import build_extension


if __name__ == "__main__":
    build_extension(force=True, verbose=True)