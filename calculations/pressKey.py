import keyboard
import time

while True:
    if keyboard.is_pressed('esc'):
        print('Bye!')
        exit()
    elif keyboard.is_pressed('q'):
        print('q')
    elif keyboard.is_pressed('r') and keyboard.is_pressed('1'):
        print('r1')
    time.sleep(0.1)
