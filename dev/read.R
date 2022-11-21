# install.packages("arrow")

library(arrow)

a <- arrow::read_feather('test.arrow')
str(a)

b <- arrow::read_feather('testc.arrow')
str(b)

c <- arrow::read_csv_arrow('test.csv')
str(c)

arrow::write_feather(a, 'testr.arrow', compression='lz4')
