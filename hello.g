Read("read.g");
RedisCommand("MULTI"); RedisCommand("LLEN hello"); RedisCommand("LPUSH hello 123"); RedisCommand("EXEC");
