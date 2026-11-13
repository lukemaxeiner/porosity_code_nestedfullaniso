CC = cc
TIFF_CFLAGS = 
TIFF_LIBS = -ltiff

CFLAGS = -O2 -Wall -Wextra -std=c11 $(TIFF_CFLAGS)
LDFLAGS = $(TIFF_LIBS) -lm

TARGET = porosity
SRC = main.c input.c tiff.c porosity_global.c porosity_field.c porosity_covariance_lag.c porosity_covariance_refwindow.c porosity_covariance_realization_prep.c porosity_covariance_tail_noise.c porosity_covariance_realization_generation.c
TEST_TARGET = porosity_tests
TEST_SRC = porosity_tests.c input.c tiff.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

$(TEST_TARGET): $(TEST_SRC) $(SRC)
	$(CC) $(CFLAGS) -o $(TEST_TARGET) $(TEST_SRC) $(LDFLAGS)

test: $(TEST_TARGET)
	./$(TEST_TARGET)

clean:
	rm -f $(TARGET) $(TEST_TARGET) *.o
