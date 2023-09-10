import os
current_path = os.path.abspath(__file__)
print(current_path)
os.chdir("//")
current_path = os.path.abspath(__file__)
print(current_path)
print(os.getcwd())
