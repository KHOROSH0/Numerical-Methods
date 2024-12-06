import time
for i in range(200):
    time.sleep(0.5*(len(str(i))))
    if i%10==1 and i != 11:
        print(f'{i} manul')
    elif (i%10 == 2 or i%10 == 3 or i%10 == 4) and (i < 10 or i > 20):
        print(f'{i} manula')
    else:
        print(f'{i} manulov')