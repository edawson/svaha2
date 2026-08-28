# Build stage
FROM debian:bookworm AS builder

# Install build dependencies
RUN apt-get update && apt-get install -y \
    clang \
    make \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /build

# Copy source code
COPY include/ include/
COPY src/ src/
COPY Makefile .

# Build the application
# We override CXXFLAGS to remove -march=native for portability in the container
RUN make CXXFLAGS="-O3 -ffast-math -std=c++17"

# Runtime stage
FROM debian:bookworm-slim

# Create a non-root user
RUN groupadd -r svaha && useradd -r -g svaha svaha

# Set working directory
WORKDIR /app

# Copy the binary from the builder stage
COPY --from=builder /build/svaha /usr/local/bin/svaha

# Set permissions
RUN chown svaha:svaha /usr/local/bin/svaha

# Switch to non-root user
USER svaha

# Set entrypoint
ENTRYPOINT ["svaha"]

# Default command
CMD ["--help"]
