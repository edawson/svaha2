# Root Makefile for svaha2 (C++ and Rust)

VERSION:=0.2.0
DOCKER_IMAGE:=erictdawson/svaha2:v$(VERSION)
DOCKER_LATEST:=erictdawson/svaha2:latest

.PHONY: all cpp rust clean test docker-build docker-push

all: cpp rust

cpp:
	$(MAKE) -C cpp
	cp cpp/svaha ./svaha-cpp

rust:
	cd svaha-rs && cargo build --release
	cp svaha-rs/target/release/svaha-rs ./svaha-rust

test: test-cpp test-rust

test-cpp:
	$(MAKE) -C cpp test

test-rust:
	cd svaha-rs && cargo test

clean:
	$(MAKE) -C cpp clean
	cd svaha-rs && cargo clean
	$(RM) svaha-cpp svaha-rust

docker-build:
	docker build -t $(DOCKER_IMAGE) -t $(DOCKER_LATEST) .

docker-push: docker-build
	docker push $(DOCKER_IMAGE)
	docker push $(DOCKER_LATEST)
