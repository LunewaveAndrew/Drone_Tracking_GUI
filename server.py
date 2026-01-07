import json
import queue
import socket
import threading
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

HOST = "0.0.0.0"
HTTP_PORT = 8000

UDP_HOST = "127.0.0.1"
UDP_PORT = 9999

clients_lock = threading.Lock()
clients = set()

sensor_state = {"lat": None, "lon": None, "t": None}


def push_msg(msg: dict):
    with clients_lock:
        for cq in list(clients):
            try:
                cq.put_nowait(msg)
            except Exception:
                pass


def udp_listener():
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.bind((UDP_HOST, UDP_PORT))

    while True:
        data, _addr = sock.recvfrom(65535)
        try:
            msg = json.loads(data.decode("utf-8").strip())

            if "az" not in msg and "final_AZ_ITP" in msg:
                msg["az"] = msg["final_AZ_ITP"]
            if "el" not in msg and "final_EL_ITP" in msg:
                msg["el"] = msg["final_EL_ITP"]

            push_msg(msg)
        except Exception:
            pass


class Handler(SimpleHTTPRequestHandler):
    def do_GET(self):
        if self.path == "/events":
            self.send_response(200)
            self.send_header("Content-Type", "text/event-stream")
            self.send_header("Cache-Control", "no-cache")
            self.end_headers()

            client_q = queue.Queue(maxsize=200)
            with clients_lock:
                clients.add(client_q)

            try:
                if sensor_state["lat"] is not None and sensor_state["lon"] is not None:
                    init = {
                        "kind": "sensor",
                        "lat": sensor_state["lat"],
                        "lon": sensor_state["lon"],
                        "t": sensor_state["t"],
                    }
                    self.wfile.write(b"data: " + json.dumps(init).encode("utf-8") + b"\n\n")
                    self.wfile.flush()

                while True:
                    msg = client_q.get()
                    payload = json.dumps(msg).encode("utf-8")
                    self.wfile.write(b"data: " + payload + b"\n\n")
                    self.wfile.flush()

            except (BrokenPipeError, ConnectionResetError):
                pass
            finally:
                with clients_lock:
                    clients.discard(client_q)
            return

        if self.path == "/sensor":
            payload = json.dumps({"ok": True, "sensor": sensor_state}).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Cache-Control", "no-cache")
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)
            return

        super().do_GET()

    def do_POST(self):
        if self.path != "/sensor":
            self.send_error(404)
            return

        try:
            length = int(self.headers.get("Content-Length", "0"))
            body = self.rfile.read(length) if length > 0 else b"{}"
            data = json.loads(body.decode("utf-8"))

            lat = float(data["lat"])
            lon = float(data["lon"])
            t = data.get("t", None)

            if not (-90.0 <= lat <= 90.0 and -180.0 <= lon <= 180.0):
                self.send_error(400)
                return

            sensor_state["lat"] = lat
            sensor_state["lon"] = lon
            sensor_state["t"] = t

            push_msg({"kind": "sensor", "lat": lat, "lon": lon, "t": t})

            payload = json.dumps({"ok": True, "sensor": sensor_state}).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Cache-Control", "no-cache")
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)
        except Exception:
            self.send_error(400)


if __name__ == "__main__":
    threading.Thread(target=udp_listener, daemon=True).start()
    ThreadingHTTPServer((HOST, HTTP_PORT), Handler).serve_forever()
