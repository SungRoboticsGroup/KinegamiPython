class MutableObject:
    def __init__(self):
        self.value = 0
    
    def __str__(self):
        return str(self.value)
    
    def increment(self):
        self.value += 1
        return self.value
    

def generator(mut : MutableObject):
    while mut.value < 10:
        yield mut.increment()
    

mut = MutableObject()
"""
print(mut)
gen = generator(mut)
print(next(gen))
print(mut)
print(mut.increment())
print(mut.increment())
print(mut.increment())
print(next(gen))"
"""

for i in generator(mut):
    print(i)
    mut.increment()