class day():
    def __init__(self,contacts, visits):
        self.contacts = contacts
        self.visits = visits
        
    def __str__(self):
        string = "Total, {:d} contacts and {:d} visits".format(self.contacts, self.visits)
        return string
        
    def __add__(self, other_list):
        contacts = self.contacts 
        visits = self.visits 
        for other in other_list:
            contacts += other.contacts
            visits += other.visits
        return day(contacts,visits)
        
        
if __name__ == "__main__":
    day1 = day(2,3)
    day2 = day(10,4)
    day3 = day(11,1)   
    day4 = day1+day2 
    day5 = day1 + [day2,day3]
    print(day1,day2,day3)
    print(day4)
    print(day5)