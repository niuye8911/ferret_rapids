import shutil, random, os

os.mkdir('queries50')
filenames = random.sample(os.listdir('./queriesnative'), 50)
for fname in filenames:
    srcpath = os.path.join('./queriesnative', fname)
    shutil.copyfile(srcpath, './queries50/'+fname)

