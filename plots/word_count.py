import re
f = open('/Users/jeffreyburggraf/PycharmProjects/TwoNeutronCorrelations/plots/past.txt')
words = 0
for line in f:
    words += len(re.findall(r'[^ ]+ +',line))

print('your thesis is {0} words'.format(words))