##########################################################
#            STAGED FIGURE PNG COMPRESSION SCRIPT
#  (c) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
##########################################################

# Shrinks staged PNGs before upload by converting them to 256-color palette PNGs
# with libimagequant (the engine behind pngquant). Files are typically 2.5-4x smaller
# with no visible change. Installed with pip (`imagequant`), so no sudo is needed.
#
# Usage (from the upload scripts, after the build):
#   python utils/compress_figures.py "$STAGE_DIR"/*.png
#
# This never fails the job: if imagequant is missing or a file can't be compressed,
# the original file is left as-is and uploaded normally.

import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

from PIL import Image


# parallel workers (each image needs well under 100 MB, so this is safe on the 10 GB VM)
MAX_WORKERS = 4



def compress_png(path):

    import imagequant

    with Image.open(path) as im:
        # already a palette PNG (e.g. compressed by an earlier job sharing this folder)
        if im.mode == "P":
            return 0, 0
        rgba = im.convert("RGBA")

    original_size = os.path.getsize(path)

    quantized = imagequant.quantize_pil_image(rgba, dithering_level=0.0, max_colors=256,
                                              min_quality=0, max_quality=100)
    # figures are fully opaque, so keep an RGB palette (no transparency chunk)
    quantized.putpalette(quantized.getpalette("RGB"), "RGB")
    quantized.info.pop("transparency", None)

    # write beside the original, then swap it in only if it's smaller
    temp_path = path + ".part"
    quantized.save(temp_path, format="PNG", optimize=True)
    new_size = os.path.getsize(temp_path)

    if new_size < original_size:
        os.replace(temp_path, path)
    else:
        os.remove(temp_path)
        new_size = original_size

    return original_size, new_size



def _safe_compress(path):
    try:
        return compress_png(path)
    except Exception as e:
        print(f"    PNG COMPRESSION SKIPPED -- {os.path.basename(path)} -- {e}")
        return 0, 0



if __name__ == "__main__":

    paths = [path for path in sys.argv[1:] if path.endswith(".png") and os.path.isfile(path)]

    if len(paths) > 0:

        try:
            import imagequant
        except ImportError:
            print("    PNG COMPRESSION SKIPPED -- imagequant not installed (pip install imagequant)")
            sys.exit(0)

        start = time.time()
        with ProcessPoolExecutor(max_workers=min(MAX_WORKERS, len(paths))) as pool:
            sizes = list(pool.map(_safe_compress, paths))

        before = sum(size[0] for size in sizes)
        after = sum(size[1] for size in sizes)
        if before > 0:
            print(f"    PNGS COMPRESSED.....{len(paths)} files, {before / 1e6:.1f} MB -> {after / 1e6:.1f} MB "
                  f"({time.time() - start:.0f} s)")
