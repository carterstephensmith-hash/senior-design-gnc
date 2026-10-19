"""Generate printable ArUco markers at exact physical sizes.

Run on a laptop with OpenCV (pip install opencv-python). Print the PNGs at
100% scale ("actual size", no fit-to-page) and measure the black square with a
ruler to confirm the size before testing.

Printing several sizes lets you sweep pixels-across-marker from a single
camera height instead of moving the camera.

Example:
    python make_markers.py --sizes-cm 8 6 4 3 2 1.5
"""
import argparse

import cv2
import numpy as np

DPI = 300


def cm_to_px(cm):
    return int(round(cm / 2.54 * DPI))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dict", default="DICT_4X4_50", help="ArUco dictionary (match the mission marker)")
    parser.add_argument("--id", type=int, default=0, help="marker ID")
    parser.add_argument("--sizes-cm", type=float, nargs="+", default=[8.0], help="black-square side lengths [cm]")
    args = parser.parse_args()

    dictionary = cv2.aruco.getPredefinedDictionary(getattr(cv2.aruco, args.dict))

    for size_cm in args.sizes_cm:
        side = cm_to_px(size_cm)
        marker = cv2.aruco.generateImageMarker(dictionary, args.id, side)

        # White quiet zone around the marker, plus a label
        margin = cm_to_px(1.0)
        page = np.full((side + 2 * margin + cm_to_px(0.8), side + 2 * margin), 255, np.uint8)
        page[margin:margin + side, margin:margin + side] = marker
        label = f"{args.dict} id {args.id}  side {size_cm:g} cm @ {DPI} dpi"
        cv2.putText(page, label, (margin, margin + side + cm_to_px(0.6)),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.8, 0, 2)

        name = f"marker_{args.dict}_id{args.id}_{size_cm:g}cm.png"
        cv2.imwrite(name, page, [cv2.IMWRITE_PNG_COMPRESSION, 9])
        print(f"wrote {name}")


if __name__ == "__main__":
    main()
