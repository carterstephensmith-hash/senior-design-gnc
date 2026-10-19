"""ArUco detection bench test for the Raspberry Pi 5.

Captures frames from the camera, runs OpenCV ArUco detection on each one, and
logs per-frame timing, whether the marker was found, and how many pixels wide
it appeared. One CSV per test case; a summary prints at the end.

Requires OpenCV >= 4.7 (ArucoDetector API). For a Pi camera module (CSI) it
uses picamera2; pass --usb for a USB webcam.

Examples:
    # Detection speed at 1080p, auto exposure, 30 s
    python aruco_bench.py --label speed_1080p --duration 30

    # Fixed 4 ms exposure, camera 2.0 m above an 8 cm marker
    python aruco_bench.py --label h200_8cm_exp4ms --exposure-us 4000 --height-m 2.0 --marker-cm 8
"""
import argparse
import csv
import os
import time

import cv2
import numpy as np


def open_camera(args):
    if args.usb:
        cap = cv2.VideoCapture(args.usb_index)
        cap.set(cv2.CAP_PROP_FRAME_WIDTH, args.width)
        cap.set(cv2.CAP_PROP_FRAME_HEIGHT, args.height)
        if args.exposure_us:
            cap.set(cv2.CAP_PROP_AUTO_EXPOSURE, 1)  # manual on V4L2; driver dependent
            cap.set(cv2.CAP_PROP_EXPOSURE, args.exposure_us / 100.0)

        def grab():
            ok, frame = cap.read()
            if not ok:
                raise RuntimeError("USB camera read failed")
            return cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

        return grab, cap.release

    from picamera2 import Picamera2

    picam2 = Picamera2()
    config = picam2.create_video_configuration(main={"size": (args.width, args.height), "format": "YUV420"})
    picam2.configure(config)
    if args.exposure_us:
        picam2.set_controls({"AeEnable": False, "ExposureTime": args.exposure_us, "AnalogueGain": args.gain})
    picam2.start()
    time.sleep(1.0)  # let exposure/white balance settle

    def grab():
        yuv = picam2.capture_array("main")
        return yuv[:args.height, :args.width]  # Y plane is the grayscale image

    return grab, picam2.stop


def marker_side_px(corners):
    pts = corners.reshape(4, 2)
    return float(np.mean(np.linalg.norm(pts - np.roll(pts, -1, axis=0), axis=1)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", required=True, help="test case name, used for output files")
    parser.add_argument("--width", type=int, default=1920)
    parser.add_argument("--height", type=int, default=1080)
    parser.add_argument("--exposure-us", type=int, default=0, help="fixed exposure [us]; 0 = auto")
    parser.add_argument("--gain", type=float, default=2.0, help="analogue gain when exposure is fixed")
    parser.add_argument("--dict", default="DICT_4X4_50", help="ArUco dictionary (match the mission marker)")
    parser.add_argument("--marker-id", type=int, default=0)
    parser.add_argument("--duration", type=float, default=20.0, help="test length [s]")
    parser.add_argument("--height-m", type=float, default=float("nan"), help="measured camera height, for the log")
    parser.add_argument("--marker-cm", type=float, default=float("nan"), help="printed marker size, for the log")
    parser.add_argument("--save-every", type=int, default=50, help="save every Nth annotated frame (0 = never)")
    parser.add_argument("--usb", action="store_true", help="use a USB webcam instead of picamera2")
    parser.add_argument("--usb-index", type=int, default=0)
    args = parser.parse_args()

    out_dir = os.path.join("results", args.label)
    os.makedirs(out_dir, exist_ok=True)

    dictionary = cv2.aruco.getPredefinedDictionary(getattr(cv2.aruco, args.dict))
    detector = cv2.aruco.ArucoDetector(dictionary, cv2.aruco.DetectorParameters())

    grab, close = open_camera(args)

    rows = []
    t_start = time.perf_counter()
    t_last_print = t_start
    frame_index = 0

    try:
        while time.perf_counter() - t_start < args.duration:
            t_frame = time.perf_counter()
            gray = grab()

            t0 = time.perf_counter()
            corners, ids, _ = detector.detectMarkers(gray)
            detect_ms = 1000 * (time.perf_counter() - t0)

            found, side_px = False, float("nan")
            if ids is not None:
                for c, marker_id in zip(corners, ids.flatten()):
                    if marker_id == args.marker_id:
                        found, side_px = True, marker_side_px(c)
                        break

            rows.append({
                "frame": frame_index,
                "t_s": t_frame - t_start,
                "detect_ms": detect_ms,
                "found": int(found),
                "side_px": side_px,
                "mean_brightness": float(gray.mean()),
            })

            if args.save_every and frame_index % args.save_every == 0:
                annotated = cv2.cvtColor(gray, cv2.COLOR_GRAY2BGR)
                if ids is not None:
                    cv2.aruco.drawDetectedMarkers(annotated, corners, ids)
                cv2.imwrite(os.path.join(out_dir, f"frame_{frame_index:05d}.jpg"), annotated)

            now = time.perf_counter()
            if now - t_last_print >= 1.0:
                recent = rows[-30:]
                print(f"{now - t_start:5.1f} s | frames {len(rows)} | "
                      f"detect {np.mean([r['detect_ms'] for r in recent]):5.1f} ms | "
                      f"found {100 * np.mean([r['found'] for r in recent]):5.1f} %")
                t_last_print = now

            frame_index += 1
    finally:
        close()

    elapsed = time.perf_counter() - t_start
    with open(os.path.join(out_dir, "frames.csv"), "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    detect = np.array([r["detect_ms"] for r in rows])
    found = np.array([r["found"] for r in rows])
    sides = np.array([r["side_px"] for r in rows if r["found"]])

    summary = {
        "label": args.label,
        "resolution": f"{args.width}x{args.height}",
        "exposure_us": args.exposure_us or "auto",
        "height_m": args.height_m,
        "marker_cm": args.marker_cm,
        "frames": len(rows),
        "loop_fps": len(rows) / elapsed,
        "detect_ms_mean": detect.mean(),
        "detect_ms_p95": np.percentile(detect, 95),
        "max_detect_fps": 1000 / detect.mean(),
        "detection_rate_pct": 100 * found.mean(),
        "side_px_mean": sides.mean() if sides.size else float("nan"),
    }

    print("\n=== Summary ===")
    for key, value in summary.items():
        print(f"{key:20s} {value:.2f}" if isinstance(value, float) else f"{key:20s} {value}")

    # One line per test case, appended, so all cases end up in one table
    summary_path = os.path.join("results", "summary.csv")
    write_header = not os.path.exists(summary_path)
    with open(summary_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=summary.keys())
        if write_header:
            writer.writeheader()
        writer.writerow(summary)


if __name__ == "__main__":
    main()
