#!/usr/bin/env python3
import os
import sys
import subprocess
import webbrowser
import time
import signal
import psutil

def kill_process_and_children(proc_pid):
    process = psutil.Process(proc_pid)
    for proc in process.children(recursive=True):
        proc.kill()
    process.kill()

def main():
    # Get the directory where the executable is located
    if getattr(sys, 'frozen', False):
        # Running as compiled executable
        app_path = os.path.dirname(sys.executable)
    else:
        # Running as script
        app_path = os.path.dirname(os.path.abspath(__file__))
    
    # Start Streamlit server
    streamlit_cmd = [sys.executable, "-m", "streamlit", "run", 
                    os.path.join(app_path, "app.py"),
                    "--server.address", "localhost",
                    "--server.port", "8501",
                    "--browser.serverAddress", "localhost",
                    "--server.headless", "true",
                    "--browser.gatherUsageStats", "false"]
    
    process = subprocess.Popen(streamlit_cmd)
    
    # Wait a moment for the server to start
    time.sleep(2)
    
    # Open web browser
    webbrowser.open("http://localhost:8501")
    
    try:
        # Keep the script running
        process.wait()
    except KeyboardInterrupt:
        # Handle clean shutdown
        kill_process_and_children(process.pid)
    finally:
        # Ensure all processes are cleaned up
        kill_process_and_children(process.pid)

if __name__ == "__main__":
    main()
