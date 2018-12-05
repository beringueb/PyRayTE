class day():
    def __init__(self,contacts, visits):
        self.contacts = contacts
        self.visits = visits
        
    def __str__(self):
        string = "Total, {:d} contacts and {:d} visits".format(self.contacts, self.visits)
        return string
        

        
        
if __name__ == "__main__":
    day1 = day(2,3)
    day2 = day(10,4)
    day3 = day(11,1)  
    day_list = [day1,day2,day3] 
    for day in day_list:
        day.visits = 3 
    for day in day_list:
        print(day)