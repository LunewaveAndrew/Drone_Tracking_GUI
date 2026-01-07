How to run with recorded data:
1. Download and unzip files
2. Open GUI_test_load_data.m in matlab, server.py and index.html in vscode
3. Open a terminal in vscode and run python server.py
4. In a browser, open localhost:8000. You should see the map, an arrow, and basic ui
5. Change path = "" in GUI_test_load_data based on your file location and run
6. Go back to the browser and you should see the arrow moving

GPS:
You can use laptop GPS from the browser, but it is innacurate. If we want better accuracy, we need to use our phones GPS. In order to do this, we need ngrok.
Why ngrok? Most browsers require a secure context (HTTPS or localhost) for Geolocation. ngrok gives you an HTTPS URL your phone can open. 

ngrok setup (one-time):
1. Go to https://ngrok.com/ and download ngrok
2. Make an account and create an authtoken 
3. Open ngrok and run .\ngrok.exe config add-authtoken <your_token>

To use:
1. Follow steps 1-4 of how to run with recorded data
2. Open ngrok and run ngrok http 8000. Look for a line that look like this: Forwarding                    https://adan-acorned-undetrimentally.ngrok-free.dev -> http://localhost:8000.
   Open the https url given on your phone. It should give you the same view as your laptop.
3. Make sure your phone's location services are enabled in settings and you should now be able use your phones GPS to set the sensor location
  - GPS (set once) will set the sensor location to your phones current GPS location
  - GPS (continous) will make the sensor location follow your phones location. Stop GPS ends this 

To run with live data the process is the exact same, just run GUI_withPython_withTCP_frame_number_Linux.m instead of GUI_test_load_data.m

Notes:
Waiting status indicates that we are waiting for packets from matlab
Tracking status is when we are continuously receiving valid packets from matlab (e.g. we are tracking the drone)
No detection status occurs when we receive an invalid signal from matlab
Disconnected status occurs when we stop receiving packets for 3.5 seconds
